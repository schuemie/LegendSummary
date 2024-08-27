library(ggplot2)
library(gridExtra)

target = 141100000 # Sitagliptin
comparator = 261100000 # Semaglutide
outcome = 3 # AMI
diagnostics <- readRDS("Diagnostics.rds")

stringToVars <- function(string) {
    parts <- as.numeric(unlist(strsplit(gsub("\\{|\\}", "", string), ",")))
    parts <- as.data.frame(matrix(parts, ncol = 5, byrow = TRUE))
    colnames(parts) <- c("min", "lower", "median", "upper", "max")
    return(parts)
}

databaseIdToFactor <- function(field, databaseIds) {
    return(factor(field, levels = (sort(databaseIds))))
}

createPerDbResultsTable <- function(target, comparator, outcome, connection) {

    databaseIds <- renderTranslateQuerySql(connection = connection,
                                              sql = "SELECT database_id FROM @schema.database WHERE is_meta_analysis = 0;",
                                              schema = schema,
                                              snakeCaseToCamelCase = TRUE)[, 1]

    # Equipoise
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
                             mutate(type = "Comparator")
    )
    vizData$databaseId <- databaseIdToFactor(vizData$databaseId, databaseIds)
    vizData$type <- factor(vizData$type, levels = sort(unique(vizData$type), decreasing = FALSE))

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5) |>
        mutate(label = if_else(is.na(minEquipoise),
                               NA,
                               sprintf("Equipoise = %0.2f", minEquipoise)))
    vizDbData$databaseId <- databaseIdToFactor(vizDbData$databaseId, databaseIds)
    labelY <- max(vizData$density) + 0.5

    plot1 <- ggplot(vizData, aes(x = preferenceScore, y = density)) +
        geom_area(aes(color = type, fill = type), position = "identity", alpha = 0.5) +
        geom_hline(yintercept = 0) +
        geom_label(aes(label = label), x = 0.5, y = labelY, vjust = 1, label.size = 0, data = vizDbData) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, labelY)) +
        scale_x_continuous("Preference score") +
        scale_fill_manual(values = c(alpha("#336B91", 0.6), alpha("#EB6622", 0.6))) +
        scale_color_manual(values = c(alpha("#336B91", 0.7), alpha("#EB6622", 0.7))) +
        facet_grid(databaseId~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    plot1


    # Covariate balance
    sql <- "
    SELECT database_id,
      PERCENTILE_DISC(ARRAY[0, 0.25,0.5,0.75,1]) WITHIN GROUP (ORDER BY std_diff_before) AS percentiles_before,
      PERCENTILE_DISC(ARRAY[0, 0.25,0.5,0.75,1]) WITHIN GROUP (ORDER BY std_diff_after) AS percentiles_after
    FROM @schema.covariate_balance
    WHERE target_id = @target
        AND comparator_id = @comparator
        AND analysis_id = 5 -- Matching
        AND outcome_id = 0
    GROUP BY database_id;
    "
    balance <- renderTranslateQuerySql(connection = connection,
                            sql = sql,
                            schema = schema,
                            target = target,
                            comparator = comparator,
                            snakeCaseToCamelCase = TRUE)
    vizData <- rbind(cbind(balance,
                           stringToVars(balance$percentilesBefore),
                           type = "Before",
                           y = 1),
                     cbind(balance,
                           stringToVars(balance$percentilesAfter),
                           type = "After matching",
                           y = .5))
    vizData$databaseId <- databaseIdToFactor(vizData$databaseId, databaseIds)
    vizData$type <- factor(vizData$type, levels = sort(unique(vizData$type), decreasing = TRUE))

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5) |>
        mutate(label = if_else(is.na(maxAbsStdDiffMean),
                                     NA,
                                     sprintf("Max ASDM = %0.2f", maxAbsStdDiffMean)))
    vizDbData$databaseId <- databaseIdToFactor(vizDbData$databaseId, databaseIds)

    boxPlotHeight <- 0.2
    balanceThreshold <- 0.1

    plot2 <- ggplot(vizData, aes(x = median, y = y)) +
        geom_rect(xmin = -balanceThreshold, ymin = 0, xmax = balanceThreshold, ymax = 1.5, alpha = 0.1, size = 0) +
        geom_vline(xintercept = 0) +
        geom_errorbarh(aes(xmin = min, xmax = max, color = type), height = boxPlotHeight*2) +
        geom_rect(aes(xmin = lower, ymin=y-boxPlotHeight, xmax = upper, ymax = y+boxPlotHeight, color = type, fill = type)) +
        geom_segment(aes(x = median, y = y-boxPlotHeight, xend = median, yend = y+boxPlotHeight, color = type), size = 1) +
        geom_label(aes(label = label), x = 0, y = 2, vjust = 1, label.size = 0, data = vizDbData) +
        coord_cartesian(xlim = c(-1, 1), ylim = c(0.25, 2)) +
        scale_x_continuous("Std. diff. means") +
        scale_color_manual(values = c(alpha("#e5ae38", 0.6), alpha("#6d904f", 0.6))) +
        scale_fill_manual(values = c(alpha("#e5ae38", 0.6), alpha("#6d904f", 0.6))) +
        facet_grid(databaseId~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    plot2


    # MDRR
    vizData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 2) |>
        mutate(label = if_else(is.na(mdrr),
                               NA,
                               sprintf("MDRR = %0.1f", mdrr)),
               type = "MDRR")
    vizData$databaseId <- databaseIdToFactor(vizDbData$databaseId, databaseIds)
    breaks <- c(0, 4, 10)
    plotMdrr <- ggplot(vizData) +
        geom_rect(xmin = 0, ymin = 0, xmax = 4, ymax = 1.5, alpha = 0.1, size = 0) +
        geom_vline(xintercept = 0) +
        geom_rect(aes(xmax = mdrr, color = type, fill = type), xmin = 0, ymin = 0.25, ymax = 1.25, alpha = 0.5) +
        geom_label(aes(label = label), x = 2, y = 2, vjust = 1, label.size = 0) +
        coord_cartesian(xlim = c(0, 10), ylim = c(0, 2)) +
        scale_x_continuous("MDRR", breaks = breaks) +
        scale_color_manual(values = "#336B91") +
        scale_fill_manual(values = "#336B91") +
        facet_grid(databaseId~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm"),
            # legend.text = element_text(color = "white")
        )
    plotMdrr


    # EASE
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
        filter(!grepl("Meta-analysis", databaseId)) |>
        mutate(type = "Negative control")
    vizData$databaseId <- databaseIdToFactor(vizData$databaseId, databaseIds)

    vizDbData <- diagnostics |>
        filter(targetId == target,
               comparatorId == comparator,
               outcomeId == outcome,
               analysisId == 5) |>
        mutate(label = sprintf("Ease = %0.2f", ease))
    vizDbData$databaseId <- databaseIdToFactor(vizDbData$databaseId, databaseIds)
    breaks <- c(0.5, 1, 2)
    plot3 <- ggplot(vizData, aes(x = logRr, y = seLogRr)) +
        geom_abline(slope = 1/1.96, linetype = "dashed") +
        geom_abline(slope = -1/1.96, linetype = "dashed") +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_point(aes(color = type), shape = 16, alpha = 0.5) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0, data = vizDbData) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(0, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "#336B91") +
        facet_grid(databaseId~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    plot3

    # Estimate
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
    vizData <- hoi |>
        filter(!grepl("Meta-analysis", databaseId)) |>
        mutate(type = "Calibrated estimate",
               label = if_else(is.na(calibratedRr), NA,
                               sprintf("%0.2f (%0.2f - %0.2f)", calibratedRr, calibratedCi95Lb, calibratedCi95Ub)))
    vizData$databaseId <- databaseIdToFactor(vizData$databaseId, databaseIds)

    breaks <- c(0.5, 1, 2)
    plot4 <- ggplot(vizData, aes(x = log(calibratedRr))) +
        geom_vline(xintercept = 0) +
        geom_point(aes(color = type), shape = 16, y = -0.5) +
        geom_errorbarh(aes(xmin = log(calibratedCi95Lb), xmax = log(calibratedCi95Ub), color = type, y = -0.5, height = 0.5)) +
        geom_label(aes(label = label), x = 0, y = 1, vjust = 1, label.size = 0) +
        coord_cartesian(xlim = log(c(0.25, 4)), ylim = c(-1, 1)) +
        scale_x_continuous("Hazard ratio", breaks = log(breaks), labels = breaks) +
        scale_color_manual(values = "black") +
        facet_grid(databaseId~., switch = "y") +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.25, "cm")
        )
    plot4


    plot <- grid.arrange(plot1, plot2, plotMdrr, plot3, plot4, ncol = 5, widths = c(1, 0.5, 0.25, 0.5, 0.5))

    ggsave("plot.png", plot = plot, width = 8, height = 5, dpi = 300)
}

