# Generate all PDFs for the symposium handouts

source("GetConnectionDetails.R")
library(dplyr)

# Connect ----------------------------------------------------------------------
connection <- connect(connectionDetails)

# Get TCOs ---------------------------------------------------------------------
sql <- "SET search_path = @schema;"
renderTranslateExecuteSql(connection, sql, schema = schema)

sql <- "
    select eoi1.exposure_name as target_name, eoi2.exposure_name as comparator_name, ooi1.outcome_name, tco1.*
    from
    (
    select target_id, comparator_id, outcome_id, count(database_id) as num_databases,
      sum(target_subjects) as sum_t,
      min(target_subjects) as min_t,
      max(target_subjects) as max_t,
      sum(comparator_subjects) as sum_c,
      min(comparator_subjects) as min_c,
      max(comparator_subjects) as max_c,
      sum(case when target_outcomes < 0 then 1 else target_outcomes end) as sum_to,
      min(case when target_outcomes < 0 then 1 else target_outcomes end) as min_to,
      max(case when target_outcomes < 0 then 1 else target_outcomes end) as max_to,
      sum(case when comparator_outcomes < 0 then 1 else comparator_outcomes end) as sum_co,
      min(case when comparator_outcomes < 0 then 1 else comparator_outcomes end) as min_co,
      max(case when comparator_outcomes < 0 then 1 else comparator_outcomes end) as max_co,
      avg(calibrated_log_rr) as avg_log_rr,
      min(calibrated_log_rr) as min_log_rr,
      max(calibrated_log_rr) as max_log_rr,
      avg(calibrated_p) as avg_p,
      min(calibrated_p) as min_p,
      max(calibrated_p) as max_p
    from cohort_method_result
    where calibrated_rr > 0
    and analysis_id = 2 /*PS matching on-treatment*/
    and outcome_id < 1000 /*not negative control*/
    group by target_id, comparator_id, outcome_id
    having count(database_id) >2  /*must have at least 3 sources with an estimate*/
    and sum(case when target_outcomes < 0 then 1 else target_outcomes end) > 20
    and sum(case when comparator_outcomes < 0 then 1 else comparator_outcomes end) > 20
    and sum(target_subjects) > 10000
    and sum(comparator_subjects) > 10000
    ) tco1
    inner join
    exposure_of_interest eoi1
    on tco1.target_id = eoi1.exposure_id
    and right(eoi1.exposure_name,4) = 'main'
    inner join
    exposure_of_interest eoi2
    on tco1.comparator_id = eoi2.exposure_id
    inner join
    outcome_of_interest ooi1
    on tco1.outcome_id = ooi1.outcome_id
    and ooi1.outcome_id not in (1, 2, 5, 10, 21, 22, 24, 35, 8, 27, 31, 32, 33, 42, 36);
"

tcos <- querySql(connection = connection,
                 sql = sql,
                 snakeCaseToCamelCase = TRUE)

# Test on a sample -----------------------------------------------------------------------
sampleSize <- 10
tcoSample <- tcos[sample.int(nrow(tcos), sampleSize), ]

oldWd <- setwd("handouts")

start <- Sys.time()
for (i in seq_len(nrow(tcoSample))) {
    tco <- tcoSample[i, ]
    fileName <- sprintf("Summary_t%d_c%d_o%d.pdf",
                        tco$targetId,
                        tco$comparatorId,
                        tco$outcomeId)

    quarto::quarto_render(input = "../Summary.qmd",
                          output_format = "pdf",
                          output_file = fileName,
                          execute_params = list(target = tco$targetId,
                                                comparator = tco$comparatorId,
                                                outcome = tco$outcomeId))
}
delta <- Sys.time() - start
message("Generating PDFs took ", signif(delta, 3), " ", attr(delta, "units"))

setwd(oldWd)

