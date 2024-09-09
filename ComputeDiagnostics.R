# Compute diagnostics for all TCOs and write to table

source("GetConnectionDetails.R")
library(dplyr)

# Connect ----------------------------------------------------------------------
connection <- connect(connectionDetails)

# Balance ----------------------------------------------------------------------
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    analysis_id,
    PERCENTILE_DISC(ARRAY[0, 0.25,0.5,0.75,1]) WITHIN GROUP (ORDER BY std_diff_before) AS percentiles_before,
    PERCENTILE_DISC(ARRAY[0, 0.25,0.5,0.75,1]) WITHIN GROUP (ORDER BY std_diff_after) AS percentiles_after,
    MAX(ABS(std_diff_after)) AS max_abs_std_diff_mean
FROM @schema.covariate_balance
WHERE outcome_id = 0
GROUP BY database_id,
    target_id,
    comparator_id,
    analysis_id;
"
balance <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   snakeCaseToCamelCase = TRUE)

stringToVars <- function(string, prefix) {
    parts <- as.numeric(unlist(strsplit(gsub("\\{|\\}", "", string), ",")))
    parts <- as.data.frame(matrix(parts, ncol = 5, byrow = TRUE))
    colnames(parts) <- paste0(prefix, c("Min", "Lower", "Median", "Upper", "Max"))
    return(parts)
}
balance <- bind_cols(
    balance,
    stringToVars(balance$percentilesBefore, "sdmBefore"),
    stringToVars(balance$percentilesAfter, "sdmAfter")
) |>
    select(-percentilesBefore, -percentilesAfter)

#saveRDS(balance, "e:/temp/LegendT2dmDiagnostics/balance.rds")
# balance <- readRDS("e:/temp/LegendT2dmDiagnostics/balance.rds")


# Equipoise --------------------------------------------------------------------
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    analysis_id,
    CASE WHEN
        target_equipoise < comparator_equipoise THEN target_equipoise
        ELSE comparator_equipoise
    END AS min_equipoise
FROM (
    SELECT database_id,
        target_id,
        comparator_id,
        analysis_id,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN target_density ELSE 0 END)/SUM(target_density) AS target_equipoise,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN comparator_density ELSE 0 END)/SUM(comparator_density) AS comparator_equipoise
    FROM @schema.preference_score_dist
    GROUP BY database_id,
        target_id,
        comparator_id,
        analysis_id
) tmp;
"

equipoise <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     min = 0.3,
                                     max = 0.7,
                                     snakeCaseToCamelCase = TRUE)
# saveRDS(equipoise, "e:/temp/LegendT2dmDiagnostics/equipoise.rds")
# equipoise <- readRDS("e:/temp/LegendT2dmDiagnostics/equipoise.rds")

# EASE -------------------------------------------------------------------------
negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                              sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                              schema = schema,
                                              snakeCaseToCamelCase = TRUE)[, 1]
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    outcome_id,
    analysis_id,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
WHERE outcome_id IN (@negative_control_ids)
    AND se_log_rr IS NOT NULL;
"
ncEstimates <- renderTranslateQuerySql(connection = connection,
                                       sql = sql,
                                       schema = schema,
                                       negative_control_ids = negativeControlIds,
                                       snakeCaseToCamelCase = TRUE)
groups <- ncEstimates |>
    filter(!grepl("Meta-analysis", databaseId)) |>
    group_by(databaseId, targetId, comparatorId, analysisId) %>%
    group_split()

computeEase <- function(group) {
    if (nrow(group) >= 5) {
        null <- EmpiricalCalibration::fitMcmcNull(group$logRr, group$seLogRr)
        ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
        row <- group |>
            head(1) |>
            select(-outcomeId, -logRr, -seLogRr) |>
            mutate(ease = !!ease$ease)
        return(row)
    } else {
        return(NULL)
    }
}
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ease <- ParallelLogger::clusterApply(cluster, groups, computeEase)
ParallelLogger::stopCluster(cluster)
ease <- bind_rows(ease)
# saveRDS(ease, "e:/temp/LegendT2dmDiagnostics/ease.rds")
# ease <- readRDS("e:/temp/LegendT2dmDiagnostics/ease.rds")

# MDRR -------------------------------------------------------------------------
alpha <- 0.05
power <- 0.8
z1MinAlpha <- qnorm(1 - alpha/2)
zBeta <- -qnorm(1 - power)
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    outcome_id,
    analysis_id,
    CAST(ABS(target_subjects) AS FLOAT) / (ABS(target_subjects) + ABS(comparator_subjects)) AS p_target,
    ABS(target_outcomes) + ABS(comparator_outcomes) AS total_outcomes
FROM @schema.cohort_method_result
WHERE ABS(target_subjects) + ABS(comparator_subjects) != 0;"
mdrr <- renderTranslateQuerySql(connection = connection,
                                sql = sql,
                                schema = schema,
                                snakeCaseToCamelCase = TRUE)
mdrr <- mdrr |>
    filter(!grepl("Meta-analysis", databaseId)) |>
    mutate(pComparator = 1 - pTarget) |>
    mutate(mdrr = exp(sqrt((zBeta + z1MinAlpha)^2/(totalOutcomes * pTarget * pComparator)))) |>
    select(databaseId, targetId, comparatorId, outcomeId, analysisId, mdrr)
# saveRDS(mdrr, "e:/temp/LegendT2dmDiagnostics/mdrr")

# Join into diagnostics tables -------------------------------------------------

# Balance only computed for analyses 5 and 6, but also apply to analyses 2 and 3,
# and 8 and 9.
balance <- bind_rows(
    balance,
    balance |>
        mutate(analysisId = analysisId - 3),
    balance |>
        mutate(analysisId = analysisId + 3)
)

# Equipoise should be the same for all, and applies to unadjusted analyses too.
# Although values appear to be the same across analyses, not all analyses have
# equal number of records, so picking one, and duplicating for all:
equipoise <- equipoise |>
    filter(analysisId == 2)
equipoise <- equipoise |>
    select(-analysisId) |>
    cross_join(tibble(analysisId = 1:9))

# Ignoring analyses with ID >= 10, since these do not seem to have all diagnostics
# computed or uploaded:
mdrr <- mdrr |>
    filter(analysisId < 10)
ease <- ease |>
    filter(analysisId < 10)


# Combine diagnostics into single table:
diagnostics <- mdrr |>
    full_join(equipoise,
               by = join_by(databaseId, targetId, comparatorId, analysisId)) |>
    full_join(ease,
              by = join_by(databaseId, targetId, comparatorId, analysisId)) |>
    full_join(balance,
              by = join_by(databaseId, targetId, comparatorId, analysisId)) |>
    filter(!grepl("Meta-analysis", databaseId)) |>
    mutate(unblind = !is.na(maxAbsStdDiffMean) &
               !is.na(minEquipoise) &
               !is.na(mdrr) &
               !is.na(ease) &
               maxAbsStdDiffMean < 0.15 &
               minEquipoise > 0.25 &
               mdrr < 4 &
               ease < 0.25)
diagnostics |>
    filter(unblind) |>
    group_by(databaseId) |>
    count()

writeLines(paste(sort(unique(diagnostics$databaseId[diagnostics$unblind])), collapse = "\", \""))

saveRDS(diagnostics, "Diagnostics.rds")

# Disconnect -------------------------------------------------------------------
disconnect(connection)

