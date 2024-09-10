# Recompute meta-analysis using diagnostics
# Assumes ComputeDiagnostics.R has been executed

source("GetConnectionDetails.R")
library(dplyr)

# Connect ----------------------------------------------------------------------
connection <- connect(connectionDetails)

# Load per-db estimates and filter by diagnostics ------------------------------
negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                              sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                              schema = schema,
                                              snakeCaseToCamelCase = TRUE)[, 1]

# Filter to analysis 2 for now:
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    outcome_id,
    analysis_id,
    target_subjects,
    comparator_subjects,
    target_days,
    comparator_days,
    target_outcomes,
    comparator_outcomes,
    rr,
    ci_95_lb,
    ci_95_ub,
    p,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
WHERE se_log_rr IS NOT NULL;
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)

diagnostics <- readRDS("Diagnostics.rds")

allResults <- estimates |>
    inner_join(diagnostics |>
                   filter(unblind) |>
                   select(targetId, comparatorId, outcomeId, analysisId, databaseId),
               by = join_by(targetId, comparatorId, outcomeId, analysisId, databaseId)) |>
    mutate(type = if_else(outcomeId %in% negativeControlIds, "Negative control", "Outcome of interest"))
diagnostics <- NULL
estimates <- NULL
negativeControlIds <- NULL
groups <- split(allResults, paste(allResults$targetId, allResults$comparatorId, allResults$analysisId), drop = TRUE)
allResults <- NULL

# group = allResults |> filter(targetId == 431100001, comparatorId == 451100001)
# group= groups[[25]]
computeGroupMetaAnalysis <- function(group,
                                     allControls) {
    # outcomeGroup = group |> filter(outcomeId == 1056)
    # outcomeGroup = outcomeGroups[[1]]
    computeSingleMetaAnalysis <- function(outcomeGroup) {
        # print(outcomeGroup$outcomeId[1])

        maRow <- outcomeGroup[1, ]
        outcomeGroup <- outcomeGroup[!is.na(outcomeGroup$seLogRr), ] # drops rows with zero events in T or C

        if (nrow(outcomeGroup) == 0) {
            maRow$targetSubjects <- 0
            maRow$comparatorSubjects <- 0
            maRow$targetDays <- 0
            maRow$comparatorDays <- 0
            maRow$targetOutcomes <- 0
            maRow$comparatorOutcomes <- 0
            maRow$rr <- NA
            maRow$ci95Lb <- NA
            maRow$ci95Ub <- NA
            maRow$p <- NA
            maRow$logRr <- NA
            maRow$seLogRr <- NA
            maRow$i2 <- NA
        } else if (nrow(outcomeGroup) == 1) {
            maRow <- outcomeGroup[1, ]
            maRow$i2 <- 0
        } else {
            maRow$targetSubjects <- sumMinCellCount(outcomeGroup$targetSubjects)
            maRow$comparatorSubjects <- sumMinCellCount(outcomeGroup$comparatorSubjects)
            maRow$targetDays <- sum(outcomeGroup$targetDays)
            maRow$comparatorDays <- sum(outcomeGroup$comparatorDays)
            maRow$targetOutcomes <- sumMinCellCount(outcomeGroup$targetOutcomes)
            maRow$comparatorOutcomes <- sumMinCellCount(outcomeGroup$comparatorOutcomes)
            meta <- meta::metagen(TE = outcomeGroup$logRr,
                                  seTE = outcomeGroup$seLogRr,
                                  sm = "RR",
                                  hakn = FALSE,
                                  control = list(stepadj=0.5, maxiter=1000))
            s <- summary(meta)
            maRow$i2 <- s$I2

            rnd <- s$random
            maRow$rr <- exp(rnd$TE)
            maRow$ci95Lb <- exp(rnd$lower)
            maRow$ci95Ub <- exp(rnd$upper)
            maRow$p <- rnd$p
            maRow$logRr <- rnd$TE
            maRow$seLogRr <- rnd$seTE
        }
        maRow$databaseId <- "Meta-analysis"
        maRow$sources <- paste(outcomeGroup$databaseId[order(outcomeGroup$databaseId)],
                               collapse = ";")
        return(maRow)
    }

    sumMinCellCount <- function(counts) {
        total <- sum(abs(counts))
        if (any(counts < 0)) {
            total <- -total
        }
        return(total)
    }

    analysisId <- group$analysisId[1]
    targetId <- group$targetId[1]
    comparatorId <- group$comparatorId[1]
    ParallelLogger::logInfo("Performing meta-analysis for target ", targetId, ", comparator ", comparatorId, ", analysis ", analysisId)
    outcomeGroups <- split(group, group$outcomeId, drop = TRUE)
    outcomeGroupResults <- lapply(outcomeGroups, computeSingleMetaAnalysis)

    groupResults <- do.call(rbind, outcomeGroupResults)

    ncs <- groupResults[groupResults$type == "Negative control", ]
    ncs <- ncs[!is.na(ncs$seLogRr), ]
    if (nrow(ncs) > 5) {
        null <- EmpiricalCalibration::fitMcmcNull(ncs$logRr, ncs$seLogRr) # calibrate CIs without synthesizing positive controls, assumes error consistent across effect sizes
        model <- EmpiricalCalibration::convertNullToErrorModel(null)
        calibratedP <- EmpiricalCalibration::calibrateP(null = null,
                                                        logRr = groupResults$logRr,
                                                        seLogRr = groupResults$seLogRr)
        calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(logRr = groupResults$logRr,
                                                                          seLogRr = groupResults$seLogRr,
                                                                          model = model)
        ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
        groupResults$calibratedP <- calibratedP$p
        groupResults$calibratedRr <- exp(calibratedCi$logRr)
        groupResults$calibratedCi95Lb <- exp(calibratedCi$logLb95Rr)
        groupResults$calibratedCi95Ub <- exp(calibratedCi$logUb95Rr)
        groupResults$calibratedLogRr <- calibratedCi$logRr
        groupResults$calibratedSeLogRr <- calibratedCi$seLogRr
        groupResults$ease <- rep(ease$ease, nrow(groupResults))
    } else {
        groupResults$calibratedP <- rep(NA, nrow(groupResults))
        groupResults$calibratedRr <- rep(NA, nrow(groupResults))
        groupResults$calibratedCi95Lb <- rep(NA, nrow(groupResults))
        groupResults$calibratedCi95Ub <- rep(NA, nrow(groupResults))
        groupResults$calibratedLogRr <- rep(NA, nrow(groupResults))
        groupResults$calibratedSeLogRr <- rep(NA, nrow(groupResults))
        groupResults$ease <- rep(NA, nrow(groupResults))
    }
    return(groupResults)
}

cluster <- ParallelLogger::makeCluster(12)
results <- ParallelLogger::clusterApply(cluster,
                                        groups,
                                        computeGroupMetaAnalysis)
ParallelLogger::stopCluster(cluster)
results <- bind_rows(results)

# Compute meta-analysis diagnostics --------------------------------------------
computeMdrrFromSe <- function(seLogRr, alpha = 0.05, power = 0.8) {
    # Based on the computation of a two-sided p-value, power can be computed as
    # power = 1-pnorm(qnorm(1 - alpha/2) - (log(mdrr) / seLogRr))/2
    # That can be translated in into:
    mdrr <- exp((qnorm(1 - alpha / 2) - qnorm(2 * (1 - power))) * seLogRr)
    return(mdrr)
}
maDiagnostics <- results |>
    mutate(mdrr = computeMdrrFromSe(seLogRr)) |>
    select(targetId,
           comparatorId,
           outcomeId,
           analysisId,
           i2,
           ease,
           mdrr) |>
    mutate(unblind = !is.na(i2) &
               !is.na(mdrr) &
               !is.na(ease) &
               i2 < 0.40 &
               mdrr < 4 &
               ease < 0.25)
saveRDS(maDiagnostics, "MaDiagnostics.rds")

# Save meta-analysis results ---------------------------------------------------
results <- results %>% mutate(databaseId = "Meta-analysis") %>%
    select(-type, -ease)
saveRDS(results, "MetaAnalysis.rds")
