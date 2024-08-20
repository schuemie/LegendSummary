# Compute diagnostics for all TCOs and write to table

source("GetConnectionDetails")
library(dplyr)

connection <- connect(connectionDetails)

# Balance ------------------------------------------------------------------------------------------

sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    analysis_id,
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


# Equipoise ----------------------------------------------------------------------------------------
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

# EASE ---------------------------------------------------------------------------------------------
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
ParallelLogger:stopCluster(cluster)
ease <- bind-rows(ease)



