source("GetConnectionDetails.R")
library(dplyr)

# Connect ----------------------------------------------------------------------
connection <- connect(connectionDetails)

# Upload diagnostics -----------------------------------------------------------
diagnostics <- readRDS("Diagnostics.rds")

# Postgres doesn't have boolean type, and cannot represent infinite:
diagnostics <- diagnostics |>
    mutate(unblind = as.integer(unblind),
           mdrr = if_else(is.infinite(mdrr), 999, mdrr))

insertTable(connection = connection,
            databaseSchema = schema,
            tableName = "diagnostics_martijn",
            data = diagnostics,
            dropTableIfExists = TRUE,
            createTable = TRUE,
            camelCaseToSnakeCase = TRUE)

maDiagnostics <- readRDS("MaDiagnostics.rds")
maDiagnostics <- maDiagnostics |>
    mutate(unblind = as.integer(unblind),
           mdrr = if_else(is.infinite(mdrr), 999, mdrr))
insertTable(connection = connection,
            databaseSchema = schema,
            tableName = "diagnostics_meta_analysis_martijn",
            data = maDiagnostics,
            dropTableIfExists = TRUE,
            createTable = TRUE,
            camelCaseToSnakeCase = TRUE)

# Upload meta-analysis ---------------------------------------------------------
maEstimates <- readRDS("MetaAnalysis.rds")
maEstimates <- maEstimates |>
    mutate(databaseId = "meta-analysis5")

# Verify we have all necessary columns:
oneRow <- renderTranslateQuerySql(connection = connection,
                                  sql = "SELECT TOP 1 * FROM @schema.cohort_method_result;",
                                  schema = schema,
                                  snakeCaseToCamelCase = TRUE)
colnames(oneRow)[!colnames(oneRow) %in% colnames(maEstimates)]
# character(0)


# Back up original cohort_method_result
sql <- "
CREATE TABLE @schema.cohort_method_result_backup_martijn AS
SELECT *
FROM @schema.cohort_method_result;
"
renderTranslateExecuteSql(connection, sql, schema = schema)

# Insert rows:
insertTable(connection = connection,
            databaseSchema = schema,
            tableName = "cohort_method_result",
            data = maEstimates,
            dropTableIfExists = FALSE,
            createTable = FALSE,
            camelCaseToSnakeCase = TRUE)

# Verify insertion:
newCount <- renderTranslateQuerySql(connection = connection,
                        sql = "SELECT COUNT(*) FROM @schema.cohort_method_result;",
                        schema = schema)[1, 1]
oldCount <- renderTranslateQuerySql(connection = connection,
                                    sql = "SELECT COUNT(*) FROM @schema.cohort_method_result_backup_martijn;",
                                    schema = schema)[1, 1]
(newCount-oldCount) == nrow(maEstimates)
# TRUE
