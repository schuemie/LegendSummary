# Generate all PDFs for the symposium handouts

source("GetConnectionDetails.R")
library(dplyr)

# Connect ----------------------------------------------------------------------
connection <- connect(connectionDetails)

# Get TCOs ---------------------------------------------------------------------
sql <- "SET search_path = @schema;"
renderTranslateExecuteSql(connection, sql, schema = schema)

sql <- "
    drop table if exists exposure_w_class;
    create temp table exposure_w_class as
    select exposure_id, substring(eoi1.exposure_name from 1 for length(eoi1.exposure_name) - 5) as ingredient_name,
      CASE WHEN eoi1.exposure_name ILIKE '%gliptin%' THEN 'DPP-4 inhibitors'
            WHEN eoi1.exposure_name ILIKE '%flozin%' THEN 'SGLT2 Inhibitors'
            WHEN eoi1.exposure_name ILIKE '%tide%' THEN 'GLP-1 Receptor Agonists'
            ELSE 'Sulfonylureas'
        END as class_name
    from exposure_of_interest eoi1
    where right(eoi1.exposure_name,4) = 'main'
    ;


    drop table if exists legend_pairs_to_review;
    create temp table legend_pairs_to_review as
    select eoi1.ingredient_name as target_ingredient_name, eoi1.class_name as target_class_name, eoi2.ingredient_name as comparator_ingredient_name, eoi2.class_name as comparator_class_name, ooi1.outcome_name, tco1.*
    from
    (
    select cm.target_id,
      cm.comparator_id,
      cm.outcome_id,
      count(cm.database_id) as num_databases,
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
    from cohort_method_result cm
    inner join diagnostics_meta_analysis_martijn dma
        on cm.target_id = dma.target_id
            and cm.comparator_id = dma.comparator_id
            and cm.outcome_id = dma.outcome_id
            and cm.analysis_id = dma.analysis_id
    where calibrated_rr > 0
    and cm.analysis_id = 2 /*PS matching on-treatment*/
    and cm.outcome_id < 1000 /*not negative control*/
    group by cm.target_id, cm.comparator_id, cm.outcome_id
    having count(cm.database_id) >2  /*must have at least 3 sources with an estimate*/
    and sum(case when target_outcomes < 0 then 1 else target_outcomes end) > 20
    and sum(case when comparator_outcomes < 0 then 1 else comparator_outcomes end) > 20
    and sum(target_subjects) > 10000
    and sum(comparator_subjects) > 10000
    ) tco1
    inner join
    exposure_w_class eoi1
    on tco1.target_id = eoi1.exposure_id
    inner join
    exposure_w_class eoi2
    on tco1.comparator_id = eoi2.exposure_id
    inner join
    outcome_of_interest ooi1
    on tco1.outcome_id = ooi1.outcome_id
    and ooi1.outcome_id not in (1, 2, 5, 10, 21, 22, 24, 35, 8, 27, 31, 32, 33, 42, 36, 7, 38, 60, 61, 29, 37, 28)
    where eoi1.class_name <> eoi2.class_name
    ;
"
executeSql(connection, sql)


tcos <- querySql(connection = connection,
                 sql = "SELECT * FROM legend_pairs_to_review",
                 snakeCaseToCamelCase = TRUE)
disconnect(connection)


# Test on a sample -----------------------------------------------------------------------
set.seed(0)
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

# Render all -------------------------------------------------------------------
oldWd <- setwd("handouts")

start <- Sys.time()
for (i in seq_len(nrow(tcos))) {
    tco <- tcos[i, ]
    fileName <- sprintf("Summary_t%d_c%d_o%d.pdf",
                        tco$targetId,
                        tco$comparatorId,
                        tco$outcomeId)
    if (!file.exists(fileName)) {
        writeLines(sprintf("*** Generating %d of %s ***", i, nrow(tcos)))
        quarto::quarto_render(input = "../Summary.qmd",
                              output_format = "pdf",
                              output_file = fileName,
                              execute_params = list(target = tco$targetId,
                                                    comparator = tco$comparatorId,
                                                    outcome = tco$outcomeId))
    }
}
delta <- Sys.time() - start
message("Generating PDFs took ", signif(delta, 3), " ", attr(delta, "units"))

setwd(oldWd)



# Merge PDFs ---------------------------------------------------------------------------------------
library(qpdf)
pdfFiles <- list.files("handouts", "*.pdf")
dir.create(file.path("handouts", "combined"))
numFiles <- length(pdfFiles)
batchSize <- ceiling(numFiles / 4.0)
for (start in seq(1, numFiles, by = batchSize)) {
    end <- min(start + batchSize - 1, numFiles)
    batch <- pdfFiles[start:end]

    pdf_combine(input = file.path("handouts", batch),
                output = file.path("handouts", "combined", sprintf("Combined%d_%d.pdf", start, end)))
}
