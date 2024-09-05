# Use some heuristics to identify missing data

source("GetConnectionDetails.R")
library(dplyr)

connection <- connect(connectionDetails)


# Everything that has an estimate must have:
# - Preference scores
# - Covariate balance
# - Negative controls

# Fetch estimates --------------------------------------------------------------
negativeControlIds <- renderTranslateQuerySql(
    connection = connection,
    sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
    schema = schema,
    snakeCaseToCamelCase = TRUE
)[, 1]

sql <- "
SELECT *
FROM @schema.cohort_method_result
WHERE se_log_rr IS NOT NULL;
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)


hoiEstimates <- estimates |>
    filter(!outcomeId %in% negativeControlIds,
           !grepl("Meta", databaseId))

tcs <- hoiEstimates |>
    distinct(databaseId, targetId, comparatorId)

# Balance ----------------------------------------------------------------------
sql <- "
SELECT DISTINCT database_id,
    target_id,
    comparator_id
FROM @schema.covariate_balance
WHERE outcome_id = 0;
"
balanceTcs <- renderTranslateQuerySql(connection = connection,
                                      sql = sql,
                                      schema = schema,
                                      snakeCaseToCamelCase = TRUE)

# mean(balanceTcs$targetId < balanceTcs$comparatorId)
# mean(tcs$targetId < tcs$comparatorId)

tcsMissingBalance <- tcs |>
    anti_join(balanceTcs, by = join_by(databaseId, targetId, comparatorId)) |>
    arrange(databaseId, targetId, comparatorId)
readr::write_csv(tcsMissingBalance, "DatabaseTargetComparatorsMissingBalance.csv")

tcsMissingBalance |>
    group_by(databaseId) |>
    count()


# Propensity models ------------------------------------------------------------
sql <- "
SELECT DISTINCT database_id,
    target_id,
    comparator_id
FROM @schema.propensity_model;
"
modelTcs <- renderTranslateQuerySql(connection = connection,
                                    sql = sql,
                                    schema = schema,
                                    snakeCaseToCamelCase = TRUE)
tcsMissingBalanceWithModel <- modelTcs |>
    anti_join(tcsMissingBalance)
tcsMissingBalanceWithModel |>
    group_by(databaseId) |>
    count()


# Preference score -------------------------------------------------------------
sql <- "
SELECT DISTINCT database_id,
    target_id,
    comparator_id
FROM @schema.preference_score_dist;
"
psTcs <- renderTranslateQuerySql(connection = connection,
                                 sql = sql,
                                 schema = schema,
                                 snakeCaseToCamelCase = TRUE)


tcsMissingPs <- tcs |>
    anti_join(psTcs, by = join_by(databaseId, targetId, comparatorId)) |>
    arrange(databaseId, targetId, comparatorId)
readr::write_csv(tcsMissingPs, "DatabaseTargetComparatorsMissingPs.csv")

tcsMissingPs |>
    group_by(databaseId) |>
    count()
