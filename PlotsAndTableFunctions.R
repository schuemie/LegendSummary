target = 141100000 # Sitagliptin
comparator = 261100000 # Semaglutide
outcome = 3 # AMI

createPerDbResultsTable <- function(target, comparator, outcome, connection) {

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

    sql <- "
    SELECT *
    FROM @schema.diagnostics
     WHERE target_id = @target
        AND comparator_id = @comparator
        AND outcome_id = @outcome
    "
    diagnostics <- renderTranslateQuerySql(connection = connection,
                            sql = sql,
                            schema = schema,
                            target = target,
                            comparator = comparator,
                            outcome = outcome,
                            snakeCaseToCamelCase = TRUE)

    sql <- "SELECT * FROM @schema.diagnostics;"

}
