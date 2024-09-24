library(ggplot2)
library(gridExtra)
library(dplyr)
library(DatabaseConnector)

# For specific subgroups, must choose specific exposure cohorts:
target = 141100000 # Sitagliptin
comparator = 261100000 # Semaglutide
outcome = 3 # AMI

target = 131100000
comparator = 141100000
outcome = 34

diagnostics <- readRDS("Diagnostics.rds")
maDiagnostics <- readRDS("MaDiagnostics.rds")
maEstimates <- readRDS("MetaAnalysis.rds")

databaseNameToFactor <- function(field, databaseName) {
    return(factor(field, levels = (sort(databaseName))))
}
renameDatabases <- function(databaseId) {
    databaseName <- case_when(databaseId == "DA_GERMANY" ~ "IQVIA DA Germany",
                              databaseId == "LPD_FRANCE" ~ "IQVIA LDP France",
                              databaseId == "OPENCLAIMS" ~ "IQVIA Open Claims",
                              databaseId == "OptumDod" ~ "Optum Clinformatics",
                              databaseId == "OptumEHR" ~ "Optum EHR",
                              databaseId == "VAOMOP" ~ "Veterans Affairs",
                              databaseId == "CCAE" ~ "Merative CCAE",
                              databaseId == "MDCD" ~ "Merative MDCD",
                              databaseId == "MDCR" ~ "Merative MDCR",
                              TRUE ~ databaseId
    )
    return(databaseName)
}
# These databases have at least one TCO-analysis pass diagnostics:
databaseIds <- c("CCAE", "DA_GERMANY", "LPD_FRANCE", "MDCD", "MDCR", "OPENCLAIMS", "OptumDod", "OptumEHR", "VAOMOP")
databaseNames <- renameDatabases(databaseIds)

createDiagnosticsTable <- function(target, comparator, outcome, connection) {

    # Estimate -------------------------------------------------------------------------------------
    sql <- "
    SELECT database_id,
      calibrated_rr,
      calibrated_ci_95_lb,
      calibrated_ci_95_ub
    FROM @schema.cohort_method_result
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 2
        AND outcome_id = @outcome;
    "
    hoi <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   target = target,
                                   comparator = comparator,
                                   outcome = outcome,
                                   snakeCaseToCamelCase = TRUE)
    databaseIdsWithValidEstimates <- hoi |>
        filter(databaseId %in% databaseIds,
               !is.na(calibratedRr)) |>
        pull(databaseId)
    vizData <- hoi |>
        filter(databaseId %in% databaseIdsWithValidEstimates) |>
        mutate(type = "Calibrated estimate",
               label = if_else(is.na(calibratedRr), NA,
                               sprintf("%0.2f (%0.2f - %0.2f)", calibratedRr, calibratedCi95Lb, calibratedCi95Ub)),
               origin = 0,
               databaseName = renameDatabases(databaseId))
    vizData <- tibble(databaseName = databaseNames) |>
        left_join(vizData, by = join_by(databaseName))
    vizData$databaseName <- databaseNameToFactor(vizData$databaseName, databaseNames)

    breaks <- c(0.5, 1, 2)
    breaksDb <- vizData |>
        filter(!is.na(calibratedRr)) |>
        select(databaseName) |>
        cross_join(tibble(x = log(breaks)))
    plotEstimate <- ggplot(vizData, aes(x = log(calibratedRr))) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = -1), data = breaksDb, color = "lightgray") +
        geom_vline(aes(xintercept = origin)) +
        geom_point(aes(color = type), shape = 16, y = -0.5) +
        geom_errorbarh(aes(xmin = log(calibratedCi95Lb), xmax = log(calibratedCi95Ub), color = type, y = -0.5, height = 0.5)) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, size = 3) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(-1, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "black", na.translate = FALSE) +
        facet_grid(databaseName~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotEstimate

    # Equipoise ------------------------------------------------------------------------------------
    sql <- "
    SELECT database_id,
      preference_score,
      target_density,
      comparator_density
    FROM @schema.preference_score_dist
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 5;
    "
    ps <- renderTranslateQuerySql(connection = connection,
                                  sql = sql,
                                  schema = schema,
                                  target = target,
                                  comparator = comparator,
                                  snakeCaseToCamelCase = TRUE)
    vizData <- bind_rows(ps |>
                             select(databaseId,
                                    preferenceScore,
                                    density = targetDensity) |>
                             mutate(type = "Target"),
                         ps |>
                             select(databaseId,
                                    preferenceScore,
                                    density = comparatorDensity) |>
                             mutate(type = "Comparator")) |>
        filter(databaseId %in% databaseIdsWithValidEstimates) |>
        mutate(databaseName = renameDatabases(databaseId))
    # Normalize density (per database)
    vizData <- vizData |> inner_join(vizData |>
                                         group_by(databaseName) |>
                                         summarise(maxDensity = max(density)),
                                     by = join_by(databaseName)) |>
        mutate(density = density / maxDensity)

    vizData <- tibble(databaseName = databaseNames) |>
        left_join(vizData, by = join_by(databaseName))
    vizData$databaseName <- databaseNameToFactor(vizData$databaseName, databaseNames)
    vizData$type <- factor(vizData$type, levels = sort(unique(vizData$type), decreasing = FALSE))

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5,
               databaseId %in% databaseIdsWithValidEstimates,
               !is.na(minEquipoise)) |>
        mutate(label = sprintf("Equipoise = %0.2f", minEquipoise),
               origin = 0,
               databaseName = renameDatabases(databaseId))
    vizDbData$databaseName <- databaseNameToFactor(vizDbData$databaseName, databaseNames)
    labelY <- max(vizData$density, na.rm = TRUE) + 1
    breaks <- c(0, 0.5, 1)
    breaksDb <- vizDbData |>
        select(databaseName) |>
        cross_join(tibble(x = breaks))

    plotEquipoise <- ggplot(vizData, aes(x = preferenceScore, y = density)) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = labelY), data = breaksDb, color = "lightgray") +
        geom_area(aes(color = type, fill = type), position = "identity", alpha = 0.5) +
        geom_hline(aes(yintercept = origin), data = vizDbData) +
        geom_label(aes(label = label), x = 0.5, y = labelY, vjust = 1, label.size = 0, size = 3, data = vizDbData) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, labelY)) +
        scale_x_continuous("Preference score", breaks = breaks) +
        scale_fill_manual(values = c(alpha("#336B91", 0.6), alpha("#EB6622", 0.6))) +
        scale_color_manual(values = c(alpha("#336B91", 0.7), alpha("#EB6622", 0.7))) +
        facet_grid(databaseName~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotEquipoise

    # Covariate balance ----------------------------------------------------------------------------
    vizData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5,
               databaseId %in% databaseIdsWithValidEstimates,
               !is.na(maxAbsStdDiffMean))
    vizData <- rbind(
        vizData |>
            select(databaseId,
                   min = asdmBeforeMin,
                   lower = asdmBeforeLower,
                   median = asdmBeforeMedian,
                   upper = asdmBeforeUpper,
                   max = asdmBeforeMax) |>
            mutate(type = "Before",
                   y = 1),
        vizData |>
            select(databaseId,
                   min = asdmAfterMin,
                   lower = asdmAfterLower,
                   median = asdmAfterMedian,
                   upper = asdmAfterUpper,
                   max = asdmAfterMax) |>
            mutate(type = "After",
                   y = 0.5)
    ) |>
        mutate(databaseName = renameDatabases(databaseId))
    vizData <- tibble(databaseName = databaseNames) |>
        left_join(vizData, join_by(databaseName))
    vizData$databaseName <- databaseNameToFactor(vizData$databaseName, databaseNames)
    vizData$type <- factor(vizData$type, levels = sort(unique(vizData$type), decreasing = TRUE))

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5,
               databaseId %in% databaseIdsWithValidEstimates,
               !is.na(maxAbsStdDiffMean)) |>
        mutate(label = sprintf("Max ASDM = %0.2f", maxAbsStdDiffMean),
               origin = 0,
               minX = 0,
               maxX = 0.15,
               median = 0,
               y = 0,
               databaseName = renameDatabases(databaseId))
    vizDbData$databaseName <- databaseNameToFactor(vizDbData$databaseName, databaseNames)

    breaks <- c(0, 0.15, 0.5, 1)
    breakLabels = c("0", "0.15", "0.5", "1")
    breaksDb <- vizDbData |>
        select(databaseName) |>
        cross_join(tibble(x = breaks))

    boxPlotHeight <- 0.2

    plotBalance <- ggplot(vizData, aes(x = median, y = y)) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = 2), data = breaksDb, color = "lightgray") +
        geom_rect(aes(xmin = minX, ymin = 0, xmax = maxX, ymax = 1.5), alpha = 0.1, size = 0, data = vizDbData) +
        geom_vline(aes(xintercept = origin), data = vizDbData) +
        geom_errorbarh(aes(xmin = min, xmax = max, color = type), height = boxPlotHeight*2) +
        geom_rect(aes(xmin = lower, ymin=y-boxPlotHeight, xmax = upper, ymax = y+boxPlotHeight, color = type, fill = type)) +
        geom_segment(aes(x = median, y = y-boxPlotHeight, xend = median, yend = y+boxPlotHeight, color = type), size = 1) +
        geom_label(aes(label = label), x = 0.5, y = 2, vjust = 1, label.size = 0, size = 3, data = vizDbData) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0.25, 2)) +
        scale_x_continuous("ASDM", breaks = breaks, labels = breakLabels) +
        scale_color_manual(values = c(alpha("#e5ae38", 0.6), alpha("#6d904f", 0.6)), na.translate = FALSE) +
        scale_fill_manual(values = c(alpha("#e5ae38", 0.6), alpha("#6d904f", 0.6)), na.translate = FALSE) +
        facet_grid(databaseName~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotBalance

    # MDRR -----------------------------------------------------------------------------------------
    vizData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 2,
               databaseId %in% databaseIdsWithValidEstimates,
               !is.na(mdrr)) |>
        mutate(label = sprintf("MDRR = %0.1f", mdrr),
               origin = 1,
               type = "MDRR",
               threshold = 4,
               databaseName = renameDatabases(databaseId))
    vizData <- vizData |>
        full_join(tibble(databaseName = databaseNames), by = join_by(databaseName))
    vizData$databaseName <- databaseNameToFactor(vizData$databaseName, databaseNames)
    breaks <- c(1, 4, 10)
    breaksDb <- vizData |>
        filter(!is.na(mdrr)) |>
        select(databaseName) |>
        cross_join(tibble(x = breaks))
    plotMdrr <- ggplot(vizData) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = 2), data = breaksDb, color = "lightgray") +
        geom_rect(aes(xmin = 1, ymin = 0, xmax = threshold, ymax = 1), alpha = 0.1, size = 0) +
        geom_vline(aes(xintercept = origin)) +
        geom_rect(aes(xmax = mdrr, color = type, fill = type), xmin = 1, ymin = 0.25, ymax = 0.75, alpha = 0.5) +
        geom_label(aes(label = label), x = 5.5, y = 2, vjust = 1, label.size = 0, size = 3) +
        coord_cartesian(xlim = c(1, 10), ylim = c(0, 2)) +
        scale_x_continuous("MDRR", breaks = breaks) +
        scale_color_manual(values = "#336B91", na.translate = FALSE) +
        scale_fill_manual(values = "#336B91", na.translate = FALSE) +
        facet_grid(databaseName~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm"),
        )
    # plotMdrr

    # EASE -----------------------------------------------------------------------------------------
    negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                                  sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                                  schema = schema,
                                                  snakeCaseToCamelCase = TRUE)[, 1]
    sql <- "
    SELECT database_id,
      log_rr,
      se_log_rr
    FROM @schema.cohort_method_result
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 2
        AND outcome_id IN (@negative_control_ids)
        AND se_log_rr IS NOT NULL;
    "
    ncs <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   target = target,
                                   comparator = comparator,
                                   negative_control_ids = negativeControlIds,
                                   snakeCaseToCamelCase = TRUE)
    vizData <- ncs |>
        filter(databaseId %in% databaseIdsWithValidEstimates) |>
        mutate(type = "Negative control",
               databaseName = renameDatabases(databaseId))
    vizData <- tibble(databaseName = databaseNames) |>
        left_join(vizData, by = join_by(databaseName))
    vizData$databaseName <- databaseNameToFactor(vizData$databaseName, databaseNames)

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5,
               databaseId %in% databaseIdsWithValidEstimates,
               !is.na(ease)) |>
        mutate(label = sprintf("EASE = %0.2f", ease),
               origin = 0,
               slope = 1/1.96,
               databaseName = renameDatabases(databaseId))
    vizDbData$databaseName <- databaseNameToFactor(vizDbData$databaseName, databaseNames)
    breaks <- c(0.5, 1, 2)
    breaksDb <- vizDbData |>
        select(databaseName) |>
        cross_join(tibble(x = log(breaks)))
    plotEase <- ggplot(vizData, aes(x = logRr, y = seLogRr)) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = 2), data = breaksDb, color = "lightgray") +
        geom_abline(aes(slope = slope, intercept = 0), linetype = "dashed", data = vizDbData) +
        geom_abline(aes(slope = -slope, intercept = 0), linetype = "dashed", data = vizDbData) +
        geom_hline(aes(yintercept = origin), data = vizDbData) +
        geom_vline(aes(xintercept = origin), data = vizDbData) +
        geom_point(aes(color = type), shape = 16, alpha = 0.5) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, size = 3, data = vizDbData) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(0, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "#336B91", na.translate = FALSE) +
        facet_grid(databaseName~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotEase

    # Combine plots --------------------------------------------------------------------------------
    plot <- grid.arrange(plotEquipoise, plotBalance, plotMdrr, plotEase, plotEstimate, ncol = 5, widths = c(1.1, 0.6, 0.38, 0.6, 0.6))
    #ggsave("plot.png", plot = plot, width = 8.5, height = 6, dpi = 300)
    return(plot)
}


createHowOftenTable <- function(target, comparator, outcome, connection) {
    sql <- "
    SELECT database_id,
      target_subjects,
      target_days,
      target_outcomes
    FROM @schema.cohort_method_result
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 2
        AND outcome_id = @outcome;
    "
    counts <- renderTranslateQuerySql(connection = connection,
                                      sql = sql,
                                      schema = schema,
                                      target = target,
                                      comparator = comparator,
                                      outcome = outcome,
                                      snakeCaseToCamelCase = TRUE)
    counts <- tibble(databaseId = databaseIds) |>
        left_join(counts, by = join_by(databaseId)) |>
        mutate(databaseName = renameDatabases(databaseId)) |>
        mutate(years = targetDays / 365.25) |>
        mutate(ir = 1000 * targetOutcomes / years) |>
        mutate(targetSubjects = format(targetSubjects, big.mark = ",", scientific = FALSE),
               years = format(round(years, 0), big.mark = ",", scientific = FALSE),
               targetOutcomes = format(round(targetOutcomes, 0), big.mark = ",", scientific = FALSE),
               ir = format(round(ir, 2), big.mark = ",", scientific = FALSE)) |>
        mutate(targetSubjects = gsub("-", "<", targetSubjects),
               years = gsub("-", "<", years),
               targetOutcomes = gsub("-", "<", targetOutcomes),
               ir = gsub("-", "<", ir)) |>
        mutate(targetSubjects = gsub("NA|^ *0|NaN", "-", targetSubjects),
               years = gsub("NA|^ *0|NaN", "-", years),
               targetOutcomes = gsub("NA|^ *0|NaN", "-", targetOutcomes),
               ir = gsub("NA|^0|NaN", "-", ir)) |>
        select(databaseName, targetSubjects, years, targetOutcomes, ir) |>
        arrange(databaseName)
    colnames(counts) <- c("Data source", "Persons exposed", "Person-time (yrs)", "Persons with outcome", "IR (/1,000 PY)")
    return(counts)
}

subgroups <- tibble(
    abbr = c("main",
             "younger-age",
             "middle-age",
             "older-age",
             "low-cvr",
             "higher-cvr",
             "female",
             "male",
             "black",
             "without-rdz",
             "with-rdz"),
    label = c("All",
              "Age < 45",
              "45 <= age < 65",
              "Age >= 65",
              "Lower cardiovascular risk",
              "Higher cardiovascular risk",
              "Female",
              "Male",
              "Black",
              "Without renal impairment",
              "With renal impairment"),
)

firstUp <- function(string) {
    substr(string, 1, 1) <- toupper(substr(string, 1, 1))
    return(string)
}

createHeader <- function(target, comparator, outcome, connection) {
    sql <- "
    SELECT exposure_name,
        CASE WHEN exposure_name ILIKE '%gliptin%' THEN 'DPP-4 inhibitors'
            WHEN exposure_name ILIKE '%flozin%' THEN 'SGLT2 Inhibitors'
            WHEN exposure_name ILIKE '%tide%' THEN 'GLP-1 Receptor Agonists'
            ELSE 'Sulfonylureas'
        END as class_name
    FROM @schema.exposure_of_interest
    WHERE exposure_id = @exposure_id;
    "
    targetNames <- renderTranslateQuerySql(connection = connection,
                                           sql = sql,
                                           schema = schema,
                                           exposure_id = target,
                                           snakeCaseToCamelCase = TRUE)
    comparatorNames <- renderTranslateQuerySql(connection = connection,
                                               sql = sql,
                                               schema = schema,
                                               exposure_id = comparator,
                                               snakeCaseToCamelCase = TRUE)
    sql <- "
    SELECT outcome_name
    FROM @schema.outcome_of_interest
    WHERE outcome_id = @outcome_id;
    "
    outcomeName <- renderTranslateQuerySql(connection = connection,
                                           sql = sql,
                                           schema = schema,
                                           outcome_id = outcome,
                                           snakeCaseToCamelCase = TRUE)[1, 1]
    # subgroup <- subgroups$label[which(sapply(subgroups$abbr, grepl, x = targetNames[1, 1], fixed = TRUE))]
    targetName <- firstUp(trimws(gsub(paste(subgroups$abbr, collapse = "|"), "", targetNames[1, 1])))
    comparatorName <- firstUp(trimws(gsub(paste(subgroups$abbr, collapse = "|"), "", comparatorNames[1, 1])))
    targetClass <- targetNames[1, 2]
    comparatorClass <- comparatorNames[1, 2]
    outcomeName <- gsub("outcome/", "", gsub("_", " ", outcomeName))
    header <- sprintf("- Target (class):\t**%s** (%s)\n- Comparator (class):\t**%s** (%s)\n- Outcome:\t**%s**",
                      targetName,
                      targetClass,
                      comparatorName,
                      comparatorClass,
                      outcomeName)
    return(header)
}

createMetaAnalysisTable <- function(target, comparator, outcome, connection) {

    diagnosticsRow <- maDiagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 2)

    # MDRR -----------------------------------------------------------------------------------------
    vizData <- diagnosticsRow |>
        mutate(label = sprintf("MDRR = %0.1f", mdrr),
               origin = 1,
               type = "MDRR",
               threshold = 4,
               databaseName = "Meta-analysis")
    plotMdrr <- ggplot(vizData) +
        geom_segment(x = 1, y = 0, xend = 1, yend = 2, color = "lightgray") +
        geom_segment(x = 4, y = 0, xend = 4, yend = 2, color = "lightgray") +
        geom_segment(x = 10, y = 0, xend = 10, yend = 2, color = "lightgray") +
        geom_rect(aes(xmin = 1, ymin = 0, xmax = threshold, ymax = 1), alpha = 0.1, size = 0) +
        geom_vline(aes(xintercept = origin)) +
        geom_rect(aes(xmax = mdrr, color = type, fill = type), xmin = 1, ymin = 0.25, ymax = 0.75, alpha = 0.5) +
        geom_label(aes(label = label), x = 5.5, y = 2, vjust = 1, label.size = 0, size = 3) +
        coord_cartesian(xlim = c(1, 10), ylim = c(0, 2)) +
        scale_x_continuous("MDRR", breaks = c(1, 4, 10)) +
        scale_color_manual(values = "#336B91", na.translate = FALSE) +
        scale_fill_manual(values = "#336B91", na.translate = FALSE) +
        facet_grid(databaseName ~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotMdrr

    # EASE -----------------------------------------------------------------------------------------
    negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                                  sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                                  schema = schema,
                                                  snakeCaseToCamelCase = TRUE)[, 1]

    vizData <- maEstimates |>
        filter(targetId == target,
               comparatorId == comparator,
               analysisId == 2,
               outcomeId %in% negativeControlIds) |>
        mutate(type = "Negative control")
    vizDbData <- diagnosticsRow |>
        mutate(label = sprintf("EASE = %0.2f", ease),
               origin = 0,
               slope = 1/1.96)
    breaks <- c(0.5, 1, 2)
    breaksDb <- vizDbData |>
        # select(databaseId) |>
        cross_join(tibble(x = log(breaks)))
    plotEase <- ggplot(vizData, aes(x = logRr, y = seLogRr)) +
        geom_segment(aes(x = x, y = 0, xend = x, yend = 0.75), data = breaksDb, color = "lightgray") +
        geom_abline(aes(slope = slope, intercept = 0), linetype = "dashed", data = vizDbData) +
        geom_abline(aes(slope = -slope, intercept = 0), linetype = "dashed", data = vizDbData) +
        geom_hline(aes(yintercept = origin), data = vizDbData) +
        geom_vline(aes(xintercept = origin), data = vizDbData) +
        geom_point(aes(color = type), shape = 16, alpha = 0.5) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, size = 3, data = vizDbData) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(0, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "#336B91", na.translate = FALSE) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotEase

    # Heterogeneity --------------------------------------------------------------------------------
    passingDiagnostics <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 2,
               databaseId %in% databaseIds) |>
        filter(mdrr < 4,
               minEquipoise > 0.25,
               ease < 0.25,
               maxAbsStdDiffMean < 0.15) |>
        pull(databaseId)
    sql <- "
    SELECT database_id,
      rr,
      ci_95_lb,
      ci_95_ub
    FROM @schema.cohort_method_result
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 2
        AND outcome_id = @outcome;
    "
    hoi <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   target = target,
                                   comparator = comparator,
                                   outcome = outcome,
                                   snakeCaseToCamelCase = TRUE)
    vizData <- hoi |>
        filter(databaseId %in% passingDiagnostics,
               !is.na(rr)) |>
        mutate(type = "Estimate",
               origin = 0,
               databaseName = renameDatabases(databaseId)) |>
        arrange(databaseName) |>
        mutate(y = -row_number() / n())

    vizDbData <- diagnosticsRow |>
        mutate(label = sprintf("I^2 == %0.2f", i2))

    plotHeterogeneity <- ggplot(vizData, aes(x = log(rr), y = y)) +
        geom_segment(x = log(0.5), y = 0, xend = log(0.5), yend = -1, color = "lightgray") +
        geom_segment(x = log(1), y = 0, xend = log(1), yend = -1, color = "lightgray") +
        geom_segment(x = log(2), y = 0, xend = log(2), yend = -1, color = "lightgray") +
        geom_vline(aes(xintercept = origin)) +
        geom_point(aes(color = type), shape = 16) +
        geom_errorbarh(aes(xmin = log(ci95Lb), xmax = log(ci95Ub), color = type, height = 0.1)) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, size = 3, parse = TRUE, data = vizDbData) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(-1, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "black", na.translate = FALSE) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotHeterogeneity

    # Estimate -------------------------------------------------------------------------------------
    vizData <- maEstimates  |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 2) |>
        mutate(type = "Calibrated estimate",
               label = if_else(is.na(calibratedRr), NA,
                               sprintf("%0.2f (%0.2f - %0.2f)", calibratedRr, calibratedCi95Lb, calibratedCi95Ub)),
               origin = 0)

    plotEstimate <- ggplot(vizData, aes(x = log(calibratedRr))) +
        geom_segment(x = log(0.5), y = 0, xend = log(0.5), yend = -1, color = "lightgray") +
        geom_segment(x = log(1), y = 0, xend = log(1), yend = -1, color = "lightgray") +
        geom_segment(x = log(2), y = 0, xend = log(2), yend = -1, color = "lightgray") +
        geom_vline(aes(xintercept = origin)) +
        geom_point(aes(color = type), shape = 16, y = -0.5) +
        geom_errorbarh(aes(xmin = log(calibratedCi95Lb), xmax = log(calibratedCi95Ub), color = type, y = -0.5, height = 0.5)) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, size = 3) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(-1, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "black", na.translate = FALSE) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    # plotEstimate

    # Combine plots --------------------------------------------------------------------------------
    plot <- grid.arrange(plotMdrr, plotEase, plotHeterogeneity, plotEstimate, ncol = 4, widths = c(1, 0.6, 0.6, 0.6))
    #ggsave("plotMa.png", plot = plot, width = 8.5, height = 1.5, dpi = 300)
    return(plot)
}
