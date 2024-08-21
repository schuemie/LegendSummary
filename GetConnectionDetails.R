library(DatabaseConnector)

# connectionDetails <- createConnectionDetails(
#     dbms = "postgresql",
#     server = paste(Sys.getenv("legendServer"), Sys.getenv("legendDatabase"), sep = "/"),
#     user = Sys.getenv("legendUser"),
#     password = Sys.getenv("legendPw")
# )
connectionDetails <- createConnectionDetails(
    dbms = "postgresql",
    server = paste(Sys.getenv("shinydbServer"), Sys.getenv("shinydbDatabase"), sep = "/"),
    user = Sys.getenv("shinydbUser"),
    password = Sys.getenv("shinydbPw")
)

schema <- "legendt2dm_drug_results"

# executeSql(connection, "COMMIT;")
