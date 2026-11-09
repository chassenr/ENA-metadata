# Assessing status of metadata on ENA for Illumina amplicon studies 

### prepare environment and read data ####

# setwd("C:/Users/chassenrueck/Documents/Additional_projects/Data_mining_commentary/Update_2025/")
setwd("/bio/Analysis_data/lost_found/tmp_hassenru/ENA_data_mining") # in conda env r-4.4.0
require(data.table)
require(reshape)
require(dplyr)
require(tidyverse)
require(janitor)
require(furrr)
require(dlfUtils)
require(R.utils)
require(UpSetR)
require(VennDiagram)
require(chisq.posthoc.test)

# save.image("R_analysis/status_ena_illumina_amplicon_update2025.Rdata")
# warning: large integers cannot be loaded from existing workspace --> save base count in Mbp instead to avoid large integers
load("R_analysis/status_ena_illumina_amplicon_update2025.Rdata")

# read ENA search results
ena.out <- fread(
  "ENA_search_tsv/amplicon_results_read_run_tsv.txt",
  h = T,
  sep = "\t",
  quote = "",
  nThread = 96
) %>% 
  clean_names() %>% 
  mutate(base_count_mbp = base_count/1000000)

# which columns are completely empty
empty.metadata <- colnames(ena.out)[apply(ena.out, 2, function(x) sum(is.na(x) | x == "") == length(x))]
length(empty.metadata)

# remove empty columns
ena.out.filt1 <- ena.out %>% 
  select(-all_of(empty.metadata))
dim(ena.out.filt1)
gc()
# 4844334 entries, 127 fields with data


### Filter ENA search result ####

# filter according to taxid
# warning: this is a dataset dependent filter --> if the analysis is repeated additional taxon id may have to be removed
remove_taxon <- c(
  "blood metagenome",
  "cellular supernatant metagenome",
  "ear metagenome",
  "eye metagenome",
  "feces metagenome",
  "gonad metagenome",
  "granuloma metagenome",
  "gut metagenome",
  "internal organ metagenome",
  "liver metagenome",
  "lung metagenome",
  "oral metagenome",
  "oral viral metagenome",
  "oral-nasopharyngeal metagenome",
  "peritoneum metagenome",
  "placenta metagenome",
  "reproductive system metagenome",
  "respiratory tract metagenome",
  "semen metagenome",
  "skin metagenome",
  "stomach metagenome",
  "upper respiratory tract metagenome",
  "urinary tract metagenome",
  "urine metagenome",
  "urogenital metagenome",
  "vaginal metagenome",
  "synthetic metagenome",
  "clinical metagenome",
  "electrolysis cell metagenome",
  "hospital metagenome",
  "medical device metagenome",
  "medicine metagenome",
  "nervous system metagenome",
  "nutrient bag metagenome",
  "Winogradsky column metagenome",
  "Mycobacterium tuberculosis;unclassified sequences"
)
ena.out.filt2 <- ena.out.filt1 %>% 
  filter(!scientific_name %in% remove_taxon) %>% 
  filter(!grepl("enrichment culture", scientific_name)) %>% 
  filter(!grepl("human", scientific_name))
table(ena.out.filt2$scientific_name)
table(grep("metagenom", ena.out.filt2$scientific_name, invert=T, ignore.case=T, value=T))
dim(ena.out.filt2)
# 3461316 runs

# write output to retrieve XREF publications and download xml (study, experiment and sample)
fwrite(
  ena.out.filt2,
  "ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt",
  quote = F,
  sep = "\t",
  row.names = F,
  nThread = 96
)

ena.assess <- ena.out.filt2
rm(ena.out.filt2)


### Add further metadata and define evaluation criteria ####

# identify DB entry point to SRA
# see: https://www.ddbj.nig.ac.jp/documents/prefix-e.html#dra
prefix.entry.db <- read.table("R_analysis/prefix_entry_DB.txt", h = F, sep = "\t")
ena.assess$entry_db <- substr(ena.assess$run_accession, 1, 3)

# list MIxS sample submissoin checklists according to https://www.ebi.ac.uk/ena/browser/checklists
mixs.checklists <- paste0("ERC0000", c(12:25, 31, 55:58))
# NCBI reporting standard is often provided as MIMS, MIGS, or similar, i.e. not appropriate for amplicons data
# however, for now I am aatributing this inconsistency to ignorance when selecting the standard 
# and not a conscious choice.
ena.assess$mixs <- grepl("MI", ena.assess$ncbi_reporting_standard, fixed = T) | ena.assess$checklist %in% mixs.checklists

# format year
ena.assess$year_created <- as.factor(gsub("-.*", "", ena.assess$first_created))

# identify GFBio
ena.assess$gfbio <- ena.assess$broker_name == "GFBIO"

# parse file size
ena.assess$file_size_mb = sapply(strsplit(ena.assess$fastq_bytes, ";"), function(x) sum(as.numeric(x)))/1000000


### Parse data from sample xml ####

plan(multicore, workers = 192)

xml_list <- list.files(path = "ENA_search_xml/sample_csv", pattern = "*.csv")
xml_scan <- future_map(
  xml_list,
  function(x) {
    tmp <- scan(
      what = "character",
      file = paste0("ENA_search_xml/sample_csv/", x),
      sep = "\n",
      quote = "",
      comment.char = ""
    )
  }
)
xml_split <- future_map(
  xml_scan,
  function(x) {
    str_split(x, pattern = "\t")
  }
)
# check that number of header entries (tag) matches number of values
xml_check <- sapply(xml_split, function(x) length(unique(sapply(x, length))))
table(xml_check) # all good (only 1 element)

# convert to 1 long table
xml_long <- future_map_dfr(
  xml_split, 
  function(x) data.frame(sample_accession = x[[2]][x[[1]] == "PRIMARY ID"], tag = x[[1]], value = x[[2]])
)

# extract all possible xml fields to identify potential information about variables of interest
sample_xml_fields <- unique(sort(xml_long$tag))

# identify fields potentially providing further information about primers, location, date
primers_xml <- grep("conc|cop|qpcr|plate|well|velocity", grep("primer|target.gene|subfragment|pcr|loc[ui]", sample_xml_fields, value = T, ignore.case = T), invert = T, value = T, ignore.case = T)
location_xy_xml <- grep("plate|well", grep("latitud|longitud|coord|[^a-zA-Z0-9]lat[^a-zA-Z0-9]|[^a-zA-Z0-9]lon[^a-zA-Z0-9]|[^a-zA-Z0-9]long[a-zA-Z0-9]|^lat$|^lon$|^long$|spatial.*location", sample_xml_fields, value = T, ignore.case = T), invert = T, value = T, ignore.case = T)
location_z_xml <- grep("tot|secchi|per|mol", grep("depth|elevation|altitude", sample_xml_fields, value = T, ignore.case = T), value = T, ignore.case = T, invert = T)
location_region_xml <- grep("country|geographic.*region|geo.loc|sampl.*loc", sample_xml_fields, value = T, ignore.case = T)
date_xml <- grep("collect.*date|collect.*time|^date$|event.*date|sampl.*da[yt]|sampl.*time|time.*point", sample_xml_fields, value = T, ignore.case = T)
# this is a rough approximation and may contain fals positives as well as false negatives
# it is however more reproducible and far better than doing it manually or asking chatGPT

# parse wide table for each set of metadata fields
df_xml_primers <- xml_long %>% 
  filter(tag %in% primers_xml) %>% 
  filter(value != "" & !is.na(value) & value != "NA")
df_xml_xy <- xml_long %>% 
  filter(tag %in% location_xy_xml) %>% 
  filter(value != "" & !is.na(value) & value != "NA")
df_xml_z <- xml_long %>% 
  filter(tag %in% location_z_xml) %>% 
  filter(value != "" & !is.na(value) & value != "NA")
df_xml_region <- xml_long %>% 
  filter(tag %in% location_region_xml) %>% 
  filter(value != "" & !is.na(value) & value != "NA")
df_xml_date <- xml_long %>% 
  filter(tag %in% date_xml) %>% 
  filter(value != "" & !is.na(value) & value != "NA")


### Parse data from experiment xml ####

plan(multicore, workers = 192)

xml_list_exp <- list.files(path = "ENA_search_xml/experiment_csv", pattern = "*.csv")
xml_scan_exp <- future_map(
  xml_list_exp,
  function(x) {
    tmp <- scan(
      what = "character",
      file = paste0("ENA_search_xml/experiment_csv/", x),
      sep = "\n",
      quote = "",
      comment.char = ""
    )
  }
)
xml_split_exp <- future_map(
  xml_scan_exp,
  function(x) {
    str_split(x, pattern = "\t")
  }
)
# check that number of header entries (tag) matches number of values
xml_check_exp <- sapply(xml_split_exp, function(x) length(unique(sapply(x, length))))
table(xml_check_exp) # all good (only 1 element)

# convert to 1 long table
xml_long_exp <- future_map_dfr(
  xml_split_exp, 
  function(x) data.frame(experiment_accession = x[[2]][x[[1]] == "PRIMARY ID"], tag = x[[1]], value = x[[2]])
)

# extract all possible xml fields to identify potential information about variables of interest
exp_xml_fields <- unique(sort(xml_long_exp$tag))

# identify fields potentially providing further information about primers
primers_xml_exp <- grep("conc|cop|qpcr|plate|well|velocity", grep("primer|target.gene|subfragment|pcr|design.description|loc[ui]", exp_xml_fields, value = T, ignore.case = T), invert = T, value = T, ignore.case = T)

# parse wide table for each set of metadata fields
df_xml_primers_exp <- xml_long_exp %>% 
  filter(tag %in% primers_xml_exp) %>% 
  filter(value != "" & !is.na(value) & value != "NA")


### basic stats ####

nrow(ena.assess)
# 3461316
length(unique(ena.assess$study_accession))
# 44421
length(unique(ena.assess$sample_accession))
# 2986834
length(unique(ena.assess$experiment_accession))
# 3455898
sum(unique(ena.assess$sample_accession) %in% gsub(".csv", "", xml_list, fixed = T)) # number of parsable XML with data
# 2986628
quantile(c(table(ena.assess$sample_accession)), seq(0.9, 1, 0.01))
# 90%  91%  92%  93%  94%  95%  96%  97%  98%  99% 100%
# 1    1    1    1    2    2    2    2    2    3 1157
quantile(c(table(ena.assess$sample_accession)), seq(0.99, 1, 0.001))
# 99% 99.1% 99.2% 99.3% 99.4% 99.5% 99.6% 99.7% 99.8% 99.9%  100%
# 3     3     4     4     4     5     6     6     9    20  1157
sum(c(table(ena.assess$sample_accession)) == 1)/length(unique(ena.assess$sample_accession))
# 93% of samples have only 1 run

hsize(sum(ena.assess$file_size_mb * 1000000))
# 117.1 TB of data
# instead/additional to counting the number of data sets for the summary plots, data volume could also be shown

head(sort(table(ena.out$broker_name), decreasing = T))
#           GFBIO  DTU-GE    CNSA    IBDC SeBiMER
# 4816375   16555    6849    1345     960     895
# GFBio is by far the most used brokerage service worldwide


### check taxid ####
sum(is.na(ena.assess$tax_id) | ena.assess$tax_id == "")
# taxid is always provided

# most frequent taxids/scientific_name
stats.taxid <- ena.assess %>% 
  group_by(scientific_name) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.taxid, "Summaries/stats_taxid.txt", quote = F, sep = "\t", row.names = F)


### check base count ####

sum(is.na(ena.assess$base_count_mbp))
sum(ena.assess$base_count_mbp == 0)
# 8855 runs without data
sum(is.na(ena.assess$file_size_mb))
sum(ena.assess$file_size_mb == 0)
# 8924 runs with 0 file size
sum(ena.assess$file_size_mb == 0 & ena.assess$base_count_mbp > 0) 
# 72 runs: could be due to technical errors, still data is not available
sum(ena.assess$file_size_mb > 0 & ena.assess$base_count_mbp == 0) 
# 3 runs --> SRR7475706 looks like genuine data, the other 2 were likely still in processing when the search was performed (first created 26.08.2025)

ena.assess$no_data <- ifelse(ena.assess$base_count_mbp == 0 | ena.assess$file_size_mb == 0, "no_data", "with_data")
ena.assess$no_data <- factor(ena.assess$no_data, levels = c("with_data", "no_data"))
stats.data_availability <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, no_data, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.data_availability, "Summaries/summary_no_data.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.data_availability, "Summaries/summary_no_data.rds")

# any bad examples from GFBio
ena.assess %>% filter(gfbio & mixs & no_data == "no_data") %>% pull(study_accession) %>% table()
# PRJEB42826 PRJEB43016 PRJEB54709  PRJEB57442 PRJEB72147
# 11         45         11          8          2


### check collection date ####

sum(is.na(ena.assess$collection_date)) # 17931
sum(ena.assess$collection_date == "", na.rm = T) # 107274

sum(is.na(ena.assess$collection_date_end)) # 0
sum(ena.assess$collection_date_end == "", na.rm = T) # 287362
sum(is.na(ena.assess$collection_date_start)) # 0
sum(ena.assess$collection_date_start == "", na.rm = T) # 287197

sum(
  (is.na(ena.assess$collection_date) | ena.assess$collection_date == "") & 
  (ena.assess$collection_date_end != "" | ena.assess$collection_date_start != "")
) # 599
# there are many cases where collection date is not provided, but a start/end time is available
# this applies to some TARA samples as well

# categories for assessment: 
#   syntax-compliant with concrete number (either collection_date or start/end), 
#   syntax-compliant without concrete number (either collection_date or start/end),
#   not syntax-compliant (both collection_date or start/end),
#   syntax-compliant with concrete number (xml),
#   syntax-compliant without concrete number (xml),
#   other entry in xml

df_xml_date$syntax_date <- grepl(
  "(^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$)",
  df_xml_date$value
)
df_xml_date$syntax_missing <- grepl(
  "(^not applicable$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$)",
  df_xml_date$value
)
tmp <- ena.assess %>% 
  select(sample_accession, collection_date, collection_date_start, collection_date_end) %>% 
  pivot_longer(cols = c(collection_date, collection_date_start, collection_date_end)) %>% 
  filter(value != "" & !is.na(value)) %>% 
  mutate(
    syntax_date = grepl(
      "(^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$)",
      value
    ),
    syntax_missing = grepl(
      "(^not applicable$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$)",
      value
    )
  )

ena.assess$collection_date_available <- factor(
  ifelse(
    ena.assess$sample_accession %in% tmp$sample_accession,
    ifelse(
      ena.assess$sample_accession %in% tmp$sample_accession[tmp$syntax_date],
      "tsv_date",
      ifelse(
        ena.assess$sample_accession %in% tmp$sample_accession[tmp$syntax_missing],
        "tsv_missing",
        "tsv_other"
      )
    ),
    ifelse(
      ena.assess$sample_accession %in% df_xml_date$sample_accession,
      ifelse(
        ena.assess$sample_accession %in% df_xml_date$sample_accession[df_xml_date$syntax_date],
        "xml_date",
        ifelse(
          ena.assess$sample_accession %in% df_xml_date$sample_accession[df_xml_date$syntax_missing],
          "xml_missing",
          "xml_other"
        )
      ),
      "no"
    )
  ),
  levels = c("no", "xml_other", "xml_missing", "xml_date", "tsv_other", "tsv_missing", "tsv_date")
)

stats.collection_date_available <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, collection_date_available, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.collection_date_available, "Summaries/summary_collection_date_available.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.collection_date_available, "Summaries/summary_collection_date_available.rds")

# how many tsv_date runs actually follow yyyy-mm-dd (no range)?
tmp <- ifelse(
  ena.assess$collection_date != "" & !is.na(ena.assess$collection_date), 
  grepl("[12][0-9]{3}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])", ena.assess$collection_date) & !grepl("/", ena.assess$collection_date, fixed = T), 
  "no value"
)
by(tmp, ena.assess$gfbio, function(x) table(x)/length(x))
# ena.assess$gfbio: FALSE
# FALSE (no day) no value   TRUE (with day)
# 0.37494967     0.03629972 0.58875060
# ena.assess$gfbio: TRUE
# FALSE (no day) no value    TRUE (with day)
# 0.042128467    0.009368685 0.948502847
tmp_check <- c(by(ena.assess$run_accession, ena.assess$sample_accession, length))
concrete_day = c(table(cut(tmp_check[unique(ena.assess$sample_accession[tmp == "TRUE"])], breaks = seq(0, 1200, 100))))
date_range = c(table(cut(tmp_check[unique(ena.assess$sample_accession[tmp == "FALSE"])], breaks = seq(0, 1200, 100))))
(date_range - concrete_day)[!(date_range == 0 & concrete_day == 0)]
# (0,100] (100,200] (200,300] (300,400] (400,500] (500,600] (600,700] (700,800] (1.1e+03,1.2e+03] 
# -638753 42        -4        -1        9         0         2         1         2 


### check lat/lon ####

sum(is.na(ena.assess$lat))
sum(is.na(ena.assess$lon))
sum(is.na(ena.assess$lon) & !is.na(ena.assess$lat))
# 925027 NAs (no cases where only 1 value is provided)
sum(ena.assess$lat == "", na.rm = T)
sum(ena.assess$lon == "", na.rm = T)
# no empty cells
sum(is.na(ena.assess$location))
sum(ena.assess$location == "", na.rm = T) # 925027
sum(is.na(ena.assess$location_start))
sum(ena.assess$location_start == "", na.rm = T) # 925027
sum(is.na(ena.assess$location_end))
sum(ena.assess$location_end == "", na.rm = T) # 925027
table(apply(ena.assess[, c("lat", "lon", "location", "location_start", "location_end")], 1, function(x) sum(is.na(x) | x == "")))
# either all present or no value provided in any of the 5 parameters
# let's continue with lat/lon

# categories for assessment: 
#   syntax-compliant with concrete number (both lat / lon), 
#   syntax-compliant without concrete number (either lat / lon),
#   not syntax-compliant (both collection_date or lat / lon),
#   syntax-compliant with concrete number (xml),
#   syntax-compliant without concrete number (xml),
#   other entry in xml

df_xml_xy$syntax_dd <- grepl(
  "(^[+-]?[0-9]+.?[0-9]{0,8}$)",
  df_xml_xy$value
)
df_xml_xy$syntax_missing <- grepl(
  "(^not applicable$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$)",
  df_xml_xy$value
)
tmp <- ena.assess %>% 
  select(sample_accession, lat, lon) %>% 
  mutate(
    syntax_dd_lat = grepl(
      "(^[+-]?[0-9]+.?[0-9]{0,8}$)",
      lat
    ),
    syntax_missing_lat = grepl(
      "(^not applicable$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$)",
      lat
    ),
    syntax_dd_lon = grepl(
      "(^[+-]?[0-9]+.?[0-9]{0,8}$)",
      lon
    ),
    syntax_missing_lon = grepl(
      "(^not applicable$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$)",
      lon
    )
  )

ena.assess$lat_lon_available <- factor(
  ifelse(
    !is.na(ena.assess$lat),
    ifelse(
      tmp$syntax_dd_lat & tmp$syntax_dd_lon,
      "tsv_dd",
      ifelse(
        tmp$syntax_missing_lat | tmp$syntax_missing_lon,
        "tsv_missing",
        "tsv_other"
      )
    ),
    ifelse(
      ena.assess$sample_accession %in% df_xml_xy$sample_accession,
      ifelse(
        ena.assess$sample_accession %in% df_xml_xy$sample_accession[df_xml_xy$syntax_dd],
        "xml_dd",
        ifelse(
          ena.assess$sample_accession %in% df_xml_xy$sample_accession[df_xml_xy$syntax_missing],
          "xml_missing",
          "xml_other"
        )
      ),
      "no"
    )
  ),
  levels = c("no", "xml_other", "xml_missing", "xml_dd", "tsv_other", "tsv_missing", "tsv_dd")
)
# never a case without decimal degree location provided
# this is unexpected since I assume that there must have been submissions with 'missing' location information
# is ENA not respecting its own syntax rules (conversion to NA)?

stats.lat_lon_available <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, lat_lon_available, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.lat_lon_available, "Summaries/summary_lat_lon_available.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.lat_lon_available, "Summaries/summary_lat_lon_available.rds")


### check country information and syntax ####

# INSDC country list
insdc_countries <- scan("R_analysis/insdc_country_list.txt", sep = "\n", what = "character")
insdc_countries_missing <- scan("R_analysis/insdc_country_list_missing.txt", sep = "\n", what = "character")

sum(is.na(ena.assess$country))
# no NAs
sum(ena.assess$country == "")
# 256856 empty cells

# categories for assessment: 
#   syntax-compliant (country provided),
#   syntax-compliant (missing),
#   not syntax compliant
#   syntax-compliant (country provided) in xml,
#   syntax-compliant (missing) in xml,
#   other entry in xml

# search for country/sea name anywhere in entry --> I am choosing this more relaxed approach here, 
# because it is allowed for the country name to be suffixed by a region specification

ena.assess$country_available <- factor(
  ifelse(
    ena.assess$country != "",
    ifelse(
      grepl(paste(insdc_countries, collapse = "|"), ena.assess$country),
      "tsv_country",
      ifelse(
        ena.assess$country %in% insdc_countries_missing,
        "tsv_missing",
        "tsv_other"
      )
    ),
    ifelse(
      ena.assess$sample_accession %in% df_xml_region$sample_accession,
      ifelse(
        ena.assess$sample_accession %in% df_xml_region$sample_accession[grepl(paste(insdc_countries, collapse = "|"), df_xml_region$value)],
        "xml_country",
        ifelse(
          ena.assess$sample_accession %in% df_xml_region$sample_accession[df_xml_region$value %in% insdc_countries_missing],
          "xml_missing",
          "xml_other"
        )
      ),
      "no"
    )
  ),
  levels = c("no", "xml_other", "xml_missing", "xml_country", "tsv_other", "tsv_missing", "tsv_country")
)

stats.country_available <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, country_available, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.country_available, "Summaries/summary_country_available.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.country_available, "Summaries/summary_country_available.rds")


### check z-coordinates ####

sum(is.na(ena.assess$altitude))
sum(is.na(ena.assess$depth))
sum(is.na(ena.assess$elevation))
# 3422644, 9214, 3368023 NAs
sum(ena.assess$altitude == "", na.rm = T) # no empty cells
sum(ena.assess$depth == "", na.rm = T) # 2712996 empty cells
sum(ena.assess$elevation == "", na.rm = T) # no empty cells
# hardly ever provided --> at least 1 value should be available

# categories for assessment: 
#   syntax-compliant with concrete number (either elevation, depth, altitude), 
#   syntax-compliant without concrete number (all of elevation, depth, altitude),
#   not syntax-compliant (all of elevation, depth, altitude),
#   syntax-compliant with concrete number (xml),
#   syntax-compliant without concrete number (xml),
#   other entry in xml

df_xml_z$syntax_number <- grepl(
  "((0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?)",
  df_xml_z$value
)
df_xml_z$syntax_missing <- grepl(
  "((^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$))",
  df_xml_z$value
)
tmp <- ena.assess %>% 
  mutate(altitude = as.character(altitude), depth = as.character(depth), elevation = as.character(elevation)) %>% 
  select(sample_accession, altitude, depth, elevation) %>% 
  pivot_longer(cols = c(altitude, depth, elevation)) %>% 
  filter(value != "" & !is.na(value)) %>% 
  mutate(
    syntax_number = grepl(
      "((0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?)",
      value
    ),
    syntax_missing = grepl(
      "((^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)|(^missing$))",
      value
    )
  )

ena.assess$z_coordinate <- factor(
  ifelse(
    ena.assess$sample_accession %in% tmp$sample_accession,
    ifelse(
      ena.assess$sample_accession %in% tmp$sample_accession[tmp$syntax_number],
      "tsv_number",
      ifelse(
        ena.assess$sample_accession %in% tmp$sample_accession[tmp$syntax_missing],
        "tsv_missing",
        "tsv_other"
      )
    ),
    ifelse(
      ena.assess$sample_accession %in% df_xml_z$sample_accession,
      ifelse(
        ena.assess$sample_accession %in% df_xml_z$sample_accession[df_xml_z$syntax_number],
        "xml_number",
        ifelse(
          ena.assess$sample_accession %in% df_xml_z$sample_accession[df_xml_z$syntax_missing],
          "xml_missing",
          "xml_other"
        )
      ),
      "no"
    )
  ),
  levels = c("no", "xml_other", "xml_missing", "xml_number", "tsv_other", "tsv_missing", "tsv_number")
)

stats.z_coordinate <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, z_coordinate, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.z_coordinate, "Summaries/summary_z_coordinate.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.z_coordinate, "Summaries/summary_z_coordinate.rds")

# how many tsv_number runs contain characters other than a number?
tmp$without_number <- gsub(" |-", "", (gsub("((0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?)", "", tmp$value)))
sum(tmp$without_number == "" & tmp$syntax_number)/sum(tmp$syntax_number)
sum(ena.assess$sample_accession %in% tmp$sample_accession[tmp$without_number == "" & tmp$syntax_number])
# 350836 / 677624 = 52% (of all tsv_number) contain only numbers
sum(ena.assess$sample_accession %in% tmp$sample_accession[tmp$without_number == "" & tmp$syntax_number] & ena.assess$gfbio)
# 15978 / 15978 (of all gfbio and tsv_number)
head(sort(table(tmp$without_number[tmp$syntax_number]), decreasing = T), 20)
sum(ena.assess$sample_accession %in% tmp$sample_accession[tmp$without_number %in% c("cm", "m", "mm", "inches", "in", "meters", "meter", "cmbf", "mbsl") & tmp$syntax_number])
# 310815 / 338405 = 92% (of all tsv_number with more than just numbers in value) are using common units

### check nominal length ####

sum(is.na(ena.assess$nominal_length))
# 3192364
sum(ena.assess$nominal_length == "", na.rm = T)
# no empty cells
quantile(ena.assess$nominal_length, seq(0, 1, 0.05), na.rm = T)
# some very weird values (incl. some larger 10kbp and less than 50bp)

# suspicious peaks at typical read lengths
ena.assess$check_nominal_length <- factor(
  ifelse(
    is.na(ena.assess$nominal_length),
    "not_available",
    ifelse(
      ena.assess$nominal_length %in% c(0, 150, 250, 300, 500, 600) | ena.assess$nominal_length > 2000 | ena.assess$nominal_length < 50,
      "suspicious",
      "provided"
    )
  ),
  levels = c("not_available", "suspicious", "provided")
)

stats.nominal_length <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, check_nominal_length, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.nominal_length, "Summaries/summary_nominal_length.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.nominal_length, "Summaries/summary_nominal_length.rds")


### check target gene ####

# target gene available at both sample and experiment level, but included in TSV output only at sample level
# this is counterintuitive since target gene is a feature of the experiment and not the sample
sum(ena.assess$target_gene != "" & ena.assess$experiment_accession %in% df_xml_primers_exp$experiment_accession[df_xml_primers_exp$tag == "target_gene"])
# at least, there is no overlap

sum(is.na(ena.assess$target_gene))
sum(is.na(ena.assess$library_construction_protocol))
# 0, 1819 cells
sum(ena.assess$target_gene == "")
sum(ena.assess$library_construction_protocol == "", na.rm = T) 
# 3377857, 3324434 empty cells
# library_construction_protocol: entries may or may not have anything to do with library construction and primers
# no consistent terminology (no controlled vocabulary used) --> this field is not machine readable

# categories for assessment: 
#   MIxS syntax compliant with pcr_primers (tsv),
#   any other entry (tsv),
#   MIxS syntax compliant with pcr_primers (xml),
#   any other entry (xml)

# MIxS syntax: "FWD:[ACGTRKSYMWBHDVN]+[; ]REV:[ACGTRKSYMWBHDVN]+"
# to be more relaxed: [ACGTRKSYMWBHDVN]{5,}

df_xml_primers$pcr_primers <- grepl(
  "[ACGTRKSYMWBHDVN]{5,}",
  df_xml_primers$value
)
df_xml_primers_exp$pcr_primers <- grepl(
  "[ACGTRKSYMWBHDVN]{5,}",
  df_xml_primers_exp$value
)
tmp <- ena.assess %>% 
  select(experiment_accession, target_gene, library_construction_protocol) %>% 
  pivot_longer(cols = c(target_gene, library_construction_protocol)) %>% 
  filter(value != "" & !is.na(value)) %>% 
  mutate(
    pcr_primers = grepl(
      "[ACGTRKSYMWBHDVN]{5,}",
      value
    )
  )

ena.assess$target_gene_info <- factor(
  ifelse(
    ena.assess$experiment_accession %in% tmp$experiment_accession[tmp$pcr_primers],
    "tsv_primers",
    ifelse(
      ena.assess$experiment_accession %in% tmp$experiment_accession,
      "tsv_other",
      ifelse(
        ena.assess$sample_accession %in% df_xml_primers$sample_accession[df_xml_primers$pcr_primers] | ena.assess$experiment_accession %in% df_xml_primers_exp$experiment_accession[df_xml_primers$pcr_primers],
        "xml_primers",
        ifelse(
          ena.assess$sample_accession %in% df_xml_primers$sample_accession | ena.assess$experiment_accession %in% df_xml_primers_exp$experiment_accession,
          "xml_other",
          "no"
        )
      )
    )
  ),
  levels = c("no", "xml_other", "xml_primers", "tsv_other", "tsv_primers")
)

stats.target_gene_info <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, target_gene_info, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.target_gene_info, "Summaries/summary_target_gene_info.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.target_gene_info, "Summaries/summary_target_gene_info.rds")

# how many runs XML entries at both leves 
sum(ena.assess$sample_accession %in% df_xml_primers$sample_accession & ena.assess$experiment_accession %in% df_xml_primers_exp$experiment_accession)
# 87241 runs of 143177 (sample level) and 3008042 (experiment level)
sum(ena.assess$sample_accession %in% df_xml_primers$sample_accession & ena.assess$experiment_accession %in% df_xml_primers_exp$experiment_accession[df_xml_primers_exp$tag != "DESIGN DESCRIPTION"])
# 3896 when disregarding DESIGN DESCRIPTION (20332 experiment level entries left)


### check inconsistencies in the metadata ####

# for those taxids which are metagenome,
# there should not be entries with library source GENOMIC
by(ena.assess$library_source, grepl("metagenome", ena.assess$scientific_name), table, useNA = "always")
# metagenome: FALSE
# GENOMIC METAGENOMIC        <NA>
#   11666        6802           0
# metagenome: TRUE
# GENOMIC METAGENOMIC        <NA>
#   730128     2712720         0

# the search was conducted for DNA and library strategy = AMPLICON
# there should not by any RT-PCR or random PCR library selection
table(ena.assess$library_selection, useNA = "always")
# PCR     RANDOM PCR     RT-PCR       <NA>
# 3401752      17494      42070          0
# even more strange, PCR was specified in the search query
# apparently ENA is not respecting its own controlled vocabulary anymore
# only a partial pattern match was performed
# this has been fixed in later versions of the ENA portal API

# also look at number of files here: there should only be 2
# https://ena-docs.readthedocs.io/en/latest/submit/fileprep/reads.html#fastq-format
table(sapply(strsplit(ena.assess$fastq_bytes, ";"), function(x) length(x)))
# 0     1      2         3
# 8924  180764 3219573   52055

ena.assess$metadata_inconsistencies <- factor(
  ifelse(
    sapply(strsplit(ena.assess$fastq_bytes, ";"), function(x) length(x)) != 2 | ena.assess$library_selection != "PCR" | (ena.assess$library_source == "GENOMIC" & grepl("metagenome", ena.assess$scientific_name)),
    "inconsistent (meta)data",
    "no inconsistencies found"
  ), 
  levels = c("inconsistent (meta)data", "no inconsistencies found")
)

stats.metadata_inconsistencies <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, metadata_inconsistencies, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.metadata_inconsistencies, "Summaries/summary_metadata_inconsistencies.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.metadata_inconsistencies, "Summaries/summary_metadata_inconsistencies.rds")


### check correct use of ontology terms ####

all.equal(ena.assess$environment_biome, ena.assess$broad_scale_environmental_context)
all.equal(ena.assess$environment_feature, ena.assess$local_environmental_context)
all.equal(ena.assess$environment_material, ena.assess$environmental_medium)
# content the same for old and new MIxS parameter names
# stick to new names

sum(is.na(ena.assess$broad_scale_environmental_context))
sum(is.na(ena.assess$local_environmental_context))
sum(is.na(ena.assess$environmental_medium))
# no NA cells
sum(ena.assess$broad_scale_environmental_context == "")
sum(ena.assess$local_environmental_context == "")
sum(ena.assess$environmental_medium == "")
# many empty cells

# download obo of latest stable ENVO release from:
# https://github.com/EnvironmentOntology/envo
# wget https://raw.githubusercontent.com/EnvironmentOntology/envo/master/envo.obo
# extract number and label
#   grep -A2 "^\[Term\]" envo.obo | sed '/^--$/d' | grep -v "^\[Term\]" > envo_terms.txt
# this is a very crude way of parsing the obo file
# warning: there is an error in the obo file that GO:0016301 does not have a name (should be: kinase activity)
#   grep -n '^$' envo_terms.txt
# GO:0016301 kinase activity
# PATO:0001739 radiation quality
#   sed -i '11426s/^$/name: kinase activity/' envo_terms.txt
#   sed -i '12392s/^$/name: radiation quality/' envo_terms.txt
envo <- data.frame(
  matrix(
    scan(what = "character", "ENVO_files/envo_terms.txt", sep = "\n"),
    ncol = 2,
    byrow = T
  )
)
colnames(envo) <- c("id", "name")
sum(grepl("id\\:", envo$id)) == nrow(envo)
sum(grepl("name\\:", envo$name)) == nrow(envo)
envo[!grepl("name\\:", envo$name), ]
# ok
envo$id2 <- gsub("^id\\: ", "", envo$id)
envo$name2 <- gsub("^name\\: ", "", envo$name)
# remove "obsolete" from name of obsolete terms
envo$name3 <- gsub("^obsolete:* ", "", envo$name2)
# convert all names to lowercase and remove non-alphanumeric characters
# this is to check for additional compatibility after exact matches
# using make_clean_names from janitor takes forever
envo$name4 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(envo$name3)))))
# modify ENVO names in ENA the same way
ena.envo <- ena.assess %>% 
  select(broad_scale_environmental_context, local_environmental_context, environmental_medium) %>% 
  mutate(
    broad2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(broad_scale_environmental_context))))),
    local2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(local_environmental_context))))),
    medium2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(environmental_medium)))))
  )
# run pattern match in bash with grep (faster)
write.table(envo, "ENVO_files/envo_formatted.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ena.envo, "ENVO_files/ena_out_envo.txt", sep = "\t", quote = F, row.names = F, col.names = F)

"""
# bash code (using the formatted ENVO names written to file above)
# matches to ENVO ID
cut -f1 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > broad_syntax.txt
cut -f2 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > local_syntax.txt
cut -f3 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > medium_syntax.txt

cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > broad_id.txt
cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > local_id.txt
cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > medium_id.txt

# matches to ENVO name
cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > broad_name.txt
cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > local_name.txt
cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > medium_name.txt

# matches to simplified ENVO name
cut -f4 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > broad_name2.txt
cut -f5 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > local_name2.txt
cut -f6 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > medium_name2.txt
"""

# read grep output
tmp.names <- matrix(
  c("ENVO_files/broad_syntax.txt", "ENVO_files/broad_id.txt", "ENVO_files/broad_name.txt", "ENVO_files/broad_name2.txt",
    "ENVO_files/local_syntax.txt", "ENVO_files/local_id.txt", "ENVO_files/local_name.txt", "ENVO_files/local_name2.txt",
    "ENVO_files/medium_syntax.txt", "ENVO_files/medium_id.txt", "ENVO_files/medium_name.txt", "ENVO_files/medium_name2.txt"),
  nrow = 3,
  ncol = 4
)
ena.in.envo <- vector("list", length = 3)
names(ena.in.envo) <- c("broad", "local", "medium")
for(i in 1:length(ena.in.envo)) {
  ena.in.envo[[i]] <- vector("list", length = 4)
  names(ena.in.envo[[i]]) <- c("syntax", "id", "name", "name2")
  for(j in 1:length(ena.in.envo[[i]])) {
    ena.in.envo[[i]][[j]] <- read.table(
      tmp.names[i, j],
      h = F,
      sep = "\t",
      stringsAsFactors = F,
      quote = "",
      comment.char = ""
    )
    ena.in.envo[[i]][[j]]$rn <- as.numeric(gsub(":.*", "", ena.in.envo[[i]][[j]][, 1]))
  }
}

# parse output
ena.assess$check_env_broad <- ifelse(ena.assess$broad_scale_environmental_context == "", "no", "filled")
ena.assess$check_env_broad[ena.in.envo$broad$name2$rn] <- "approx"
ena.assess$check_env_broad[unique(c(ena.in.envo$broad$name$rn, ena.in.envo$broad$id$rn))] <- "name_id"
ena.assess$check_env_broad[ena.in.envo$broad$syntax$rn] <- "syntax"
ena.assess$check_env_broad <- factor(ena.assess$check_env_broad, levels = rev(c("syntax", "name_id", "approx", "filled", "no")))

ena.assess$check_env_local <- ifelse(ena.assess$local_environmental_context == "", "no", "filled")
ena.assess$check_env_local[ena.in.envo$local$name2$rn] <- "approx"
ena.assess$check_env_local[unique(c(ena.in.envo$local$name$rn, ena.in.envo$local$id$rn))] <- "name_id"
ena.assess$check_env_local[ena.in.envo$local$syntax$rn] <- "syntax"
ena.assess$check_env_local <- factor(ena.assess$check_env_local, levels = rev(c("syntax", "name_id", "approx", "filled", "no")))

ena.assess$check_env_medium <- ifelse(ena.assess$environmental_medium == "", "no", "filled")
ena.assess$check_env_medium[ena.in.envo$medium$name2$rn] <- "approx"
ena.assess$check_env_medium[unique(c(ena.in.envo$medium$name$rn, ena.in.envo$medium$id$rn))] <- "name_id"
ena.assess$check_env_medium[ena.in.envo$medium$syntax$rn] <- "syntax"
ena.assess$check_env_medium <- factor(ena.assess$check_env_medium, levels = rev(c("syntax", "name_id", "approx", "filled", "no")))

stats.check_env_broad <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, check_env_broad, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.check_env_broad, "Summaries/summary_check_env_broad.txt", sep = "\t", quote = F, row.names = F)

stats.check_env_local <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, check_env_local, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.check_env_local, "Summaries/summary_check_env_local.txt", sep = "\t", quote = F, row.names = F)

stats.check_env_medium <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, check_env_medium, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.check_env_medium, "Summaries/summary_check_env_medium.txt", sep = "\t", quote = F, row.names = F)

saveRDS(stats.check_env_broad, "Summaries/summary_check_env_broad.rds")
saveRDS(stats.check_env_local, "Summaries/summary_check_env_local.rds")
saveRDS(stats.check_env_medium, "Summaries/summary_check_env_medium.rds")


### check status environmental sample ####

table(ena.assess$environmental_sample, useNA = "ifany")
# FALSE    TRUE    <NA>
#  4070   29696 3427550
# before the ENA API update in 2023, there were no NA calls
# now, there are 3427550 NA cells (99%)
# I assume that the entries for this parameter were revised by ENA
# at the moment, parameter is as useless for data mining as before (although less false information included)

# https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/tag_querying.html
# https://academic.oup.com/nar/article/53/D1/D49/7903379
# according to ENA, this is the metadata field that should be used to identify env origin
# it is a derived value based either on taxonomy or coordinates
table(ena.assess$tag, useNA = "ifany")
table(grepl("env_", ena.assess$tag))
by(grepl("env_", ena.assess$tag), ena.assess$environmental_sample, table)
# not consistent

# confidence scheme for env tag:
assign.tag.conf <- function(x) {
  if(length(x) > 0) {
    if(nrow(x) == 1) {
      conf = "medium" # if only 1 value available, there cannot be both tax and geo tags, i.e. medium
    } else {
      if(length(unique(x$X1)) == 1) {
        conf = "low" # if there is more than 1 value, but only from 1 of tax or geo, low
      } else {
        if(length(unique(x$X2)) == 1) {
          conf = "high" # if there is both tax and geo, but only 1 env, high
        } else {
          if(nrow(x) == 2) {
            conf = "conflict" # if there is both tax and geo, and more than 1 unique env, but only 2 values, it has to be a conflict
          } else {
            a = c(by(x$X2, x$X1, c))
            b = intersect(a[[1]], a[[2]]) 
            if(length(b) == 0) {
              conf = "conflict" # if there is both tax and geo and they have nothing in common, it is a conflict
            } else {
              if(min(table(x$X1)) == 1) {
                conf = "high" # if there is both tax and geo, one of which is only represented by 1 value, and they have something in common, it is high
              } else {
                conf = "medium" # if there is both tax and geo, both of which are represented by more than 1 value, and they have something in common, it is medium
              }
            }
          }
        }
      }
    }
  } else {
    conf = "no env tag"
  }
  return(conf)
}

tag.split <- lapply(strsplit(ena.assess$tag, ";"), function(x) grep(":", x, value = T, fixed = T))
tag.df <- future_map(
  tag.split,
  function(x) {
    if(length(x) > 0) {
      data.frame(do.call("rbind", strsplit(x, ":")))
    } else {
      x
    }
  }
)
tag.conf <- future_map_chr(tag.df, assign.tag.conf)

ena.assess$env_tag_confidence <- factor(tag.conf, levels = c("no env tag", "conflict", "low", "medium", "high"))

stats.env_tag_confidence <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, env_tag_confidence, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.env_tag_confidence, "Summaries/summary_env_tag_confidence.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.env_tag_confidence, "Summaries/summary_env_tag_confidence.rds")


### check publication available ####

# studies with publication
xref.out <- unique(sort(unlist(read.table("ENA_xref_out/xref_results_primary_study_pmc.txt", sep = "\t", h = F)[, 6:7])))[-1]
pmc.out <- scan(what = "character", file = "R_analysis/pmc_matches_study_accession.txt", sep = "\n")
s2orc.out <- scan(what = "character", file = "R_analysis/s2orc_matches_study_accession.txt", sep = "\n")

tmp.study <- data.frame(unique(ena.assess[, c("study_accession", "secondary_study_accession")]))
tmp <- list(
  xref = tmp.study$study_accession[tmp.study$study_accession %in% xref.out | tmp.study$secondary_study_accession %in% xref.out], 
  pmc = tmp.study$study_accession[tmp.study$study_accession %in% pmc.out | tmp.study$secondary_study_accession %in% pmc.out],
  s2orc = tmp.study$study_accession[tmp.study$study_accession %in% s2orc.out | tmp.study$secondary_study_accession %in% s2orc.out]
)
venn.diagram(tmp, "Summaries/publication_study_level.png", imagetype = "png")

ena.assess$publication_available <- factor(
  ifelse(
    ena.assess$study_accession %in% xref.out | ena.assess$secondary_study_accession %in% xref.out, 
    "xref_found", 
    ifelse(
      ena.assess$study_accession %in% pmc.out | ena.assess$secondary_study_accession %in% pmc.out,
      "pmc_found",
      ifelse(
        ena.assess$study_accession %in% s2orc.out | ena.assess$secondary_study_accession %in% s2orc.out,
        "s2orc_found",
        "not_found" 
      )
    )
  ),
  levels = c("not_found", "s2orc_found", "pmc_found", "xref_found")
)

stats.publication_available <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, publication_available, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.publication_available, "Summaries/summary_publication_available.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.publication_available, "Summaries/summary_publication_available.rds")


### overall quality assessment ####

# define overall quality in 3 categories
machine_readable <- ena.assess$no_data == "with_data" &
  ena.assess$collection_date_available == "tsv_date" & 
  ena.assess$lat_lon_available == "tsv_dd" &
  ena.assess$metadata_inconsistencies == "no inconsistencies found" &
  (ena.assess$env_tag_confidence == "high" | (ena.assess$check_env_broad == "syntax" & ena.assess$check_env_local == "syntax" & ena.assess$check_env_medium == "syntax")) &
  ena.assess$target_gene_info %in% c("tsv_primers", "tsv_other")
   
missing_inconsistent <- (ena.assess$no_data == "no_data" |
  ena.assess$collection_date_available == "no" |
  ena.assess$lat_lon_available == "no" |
  ena.assess$target_gene_info == "no" |
  ena.assess$metadata_inconsistencies == "inconsistent (meta)data" |
  (ena.assess$env_tag_confidence == "no env tag" & (ena.assess$check_env_broad == "no" | ena.assess$check_env_local == "no" | ena.assess$check_env_medium == "no"))) &
  ena.assess$publication_available == "not_found"

ena.assess$overall_quality <- factor(
  ifelse(
    missing_inconsistent,
    "missing_inconsistent",
    ifelse(
      machine_readable,
      "machine_readable",
      "human_readable"
    )
  ), 
  levels = c("missing_inconsistent", "human_readable", "machine_readable")
)

stats.overall_quality <- ena.assess %>% 
  group_by(entry_db, mixs, gfbio, overall_quality, year_created, .drop = F) %>% 
  summarize(
    n_data = n(),
    mbp_data = sum(file_size_mb)
  )
write.table(stats.overall_quality, "Summaries/summary_overall_quality.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.overall_quality, "Summaries/summary_overall_quality.rds")


### chi-square statistic ####

# general trends: use of MIxS between entry DB
chisq.general <- chisq.test(ena.assess$entry_db, ena.assess$mixs)
chisq.general.posthoc <- chisq.posthoc.test(chisq.general$observed, method = "bonferroni")
# all significant
# DDBJ and ENA larger number of MIxS submission

# Effect of MIxS, entry DB, GFBio on metadata quality
test.vars <- c(
  "collection_date_available",
  "lat_lon_available",
  "country_available",
  "z_coordinate",
  "metadata_inconsistencies",
  "check_env_broad",
  "check_env_local",
  "check_env_medium",
  "env_tag_confidence",
  "publication_available",
  "target_gene_info"
)

# first, run the default test for all to identify those with warning
# then, group categories to avoid warnings if applicable

chisq.entry.db <- lapply(test.vars, function(x) chisq.test(ena.assess$entry_db, ena.assess %>% pull(x)))
names(chisq.entry.db) <- test.vars
sapply(chisq.entry.db, function(x) any(x$expected < 5))
# only collection date requires merging of 2 categories: xml_missing and xml_date will be combined to xml_syntax
tmp <- ena.assess %>% select(entry_db, collection_date_available)
levels(tmp$collection_date_available) <- c("no", "xml_other", "xml_syntax", "xml_syntax", "tsv_other", "tsv_missing", "tsv_date")
chisq.entry.db$collection_date_available <- chisq.test(tmp$entry_db, tmp$collection_date_available)

chisq.mixs <- lapply(test.vars, function(x) chisq.test(ena.assess$mixs, ena.assess %>% pull(x)))
names(chisq.mixs) <- test.vars
sapply(chisq.mixs, function(x) any(x$expected < 5))

chisq.gfbio <- lapply(test.vars, function(x) chisq.test(ena.assess$gfbio, ena.assess %>% pull(x)))
names(chisq.gfbio) <- test.vars
sapply(chisq.gfbio, function(x) any(x$expected < 5))
# only collection date requires merging of 2 categories: xml_missing and xml_date will be combined to xml_syntax
tmp <- ena.assess %>% select(gfbio, collection_date_available)
levels(tmp$collection_date_available) <- c("no", "xml_other", "xml_syntax", "xml_syntax", "tsv_other", "tsv_missing", "tsv_date")
chisq.gfbio$collection_date_available <- chisq.test(tmp$gfbio, tmp$collection_date_available)

# parse results and apply p-value correction
chisq.entry.db.main <- map_dfr(
  chisq.entry.db,
  function(x) {
    tibble(chisq = x$statistic, df = x$parameter, p = x$p.value)
  }
) %>% 
  mutate(metadata_field = test.vars, .before = 1) %>% 
  mutate(padj = p.adjust(p, "bonferroni"))

chisq.mixs.main <- map_dfr(
  chisq.mixs,
  function(x) {
    tibble(chisq = x$statistic, df = x$parameter, p = x$p.value)
  }
) %>% 
  mutate(metadata_field = test.vars, .before = 1) %>% 
  mutate(padj = p.adjust(p, "bonferroni"))

chisq.gfbio.main <- map_dfr(
  chisq.gfbio,
  function(x) {
    tibble(chisq = x$statistic, df = x$parameter, p = x$p.value)
  }
) %>% 
  mutate(metadata_field = test.vars, .before = 1) %>% 
  mutate(padj = p.adjust(p, "bonferroni"))

# posthoc tests
# I decided to p-value correct across all posthocs per factor

chisq.entry.db.posthoc <- lapply(chisq.entry.db, function(x) chisq.posthoc.test(x$observed, method = "none"))
chisq.mixs.posthoc <- lapply(chisq.mixs, function(x) chisq.posthoc.test(x$observed, method = "none"))
chisq.gfbio.posthoc <- lapply(chisq.gfbio, function(x) chisq.posthoc.test(x$observed, method = "none"))

chisq.entry.db.posthoc.summary <- map_dfr(
  names(chisq.entry.db.posthoc),
  function(x) {
    chisq.entry.db.posthoc[[x]] %>% 
      mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
      column_to_rownames("rn") %>%
      select(-1, -2) %>% t() %>%
      gsub("*", "", ., fixed = T) %>% 
      as.data.frame() %>% 
      rownames_to_column(., var = "categories") %>% 
      mutate(metadata_field = x, .before = 1)
  }
) %>% 
  pivot_longer(cols = -c(1:2)) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(metadata_field, categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals),
    p_adj = p.adjust(p_values, "bonferroni")
  )

chisq.mixs.posthoc.summary <- map_dfr(
  names(chisq.mixs.posthoc),
  function(x) {
    chisq.mixs.posthoc[[x]] %>% 
      mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
      column_to_rownames("rn") %>%
      select(-1, -2) %>% t() %>%
      gsub("*", "", ., fixed = T) %>% 
      as.data.frame() %>% 
      rownames_to_column(., var = "categories") %>% 
      mutate(metadata_field = x, .before = 1)
  }
) %>% 
  pivot_longer(cols = -c(1:2)) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(metadata_field, categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals),
    p_adj = p.adjust(p_values, "bonferroni")
  )

chisq.gfbio.posthoc.summary <- map_dfr(
  names(chisq.gfbio.posthoc),
  function(x) {
    chisq.gfbio.posthoc[[x]] %>% 
      mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
      column_to_rownames("rn") %>%
      select(-1, -2) %>% t() %>%
      gsub("*", "", ., fixed = T) %>% 
      as.data.frame() %>% 
      rownames_to_column(., var = "categories") %>% 
      mutate(metadata_field = x, .before = 1)
  }
) %>% 
  pivot_longer(cols = -c(1:2)) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(metadata_field, categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals),
    p_adj = p.adjust(p_values, "bonferroni")
  )

# overall quality assessment

chisq.entry.db.oq <- chisq.test(ena.assess$entry_db, ena.assess$overall_quality)
chisq.entry.db.oq.posthoc <- chisq.posthoc.test(chisq.entry.db.oq$observed, method = "bonferroni") %>% 
  mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
  column_to_rownames("rn") %>%
  select(-1, -2) %>% t() %>%
  gsub("*", "", ., fixed = T) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "categories") %>% 
  pivot_longer(cols = -1) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals)
  )

chisq.mixs.oq <- chisq.test(ena.assess$mixs, ena.assess$overall_quality)
chisq.mixs.oq.posthoc <- chisq.posthoc.test(chisq.mixs.oq$observed, method = "bonferroni") %>% 
  mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
  column_to_rownames("rn") %>%
  select(-1, -2) %>% t() %>%
  gsub("*", "", ., fixed = T) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "categories") %>% 
  pivot_longer(cols = -1) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals)
  )

chisq.gfbio.oq <- chisq.test(ena.assess$gfbio, ena.assess$overall_quality)
chisq.gfbio.oq.posthoc <- chisq.posthoc.test(chisq.gfbio.oq$observed, method = "bonferroni") %>% 
  mutate(rn = paste(Dimension, Value, sep = "___")) %>% 
  column_to_rownames("rn") %>%
  select(-1, -2) %>% t() %>%
  gsub("*", "", ., fixed = T) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "categories") %>% 
  pivot_longer(cols = -1) %>% separate(name, into = c("group", "statistic"), sep = "___") %>% 
  pivot_wider(id_cols = c(categories, group), names_from = statistic, values_from = value) %>% 
  clean_names() %>% 
  mutate(
    p_values = as.numeric(p_values), 
    residuals = as.numeric(residuals)
  )

# write output

write.table(
  bind_rows(
    chisq.entry.db.main %>% mutate(factor = "entry DB", .before = 1),
    chisq.mixs.main %>% mutate(factor = "MIxS", .before = 1),
    chisq.gfbio.main %>% mutate(factor = "GFBio", .before = 1)
  ),
  "Summaries/chisq_main.txt",
  quote = F,
  sep = "\t",
  row.names = F
)
write.table(
  bind_rows(
    chisq.entry.db.posthoc.summary %>% mutate(factor = "entry DB", .before = 1),
    chisq.mixs.posthoc.summary %>% mutate(factor = "MIxS", .before = 1),
    chisq.gfbio.posthoc.summary %>% mutate(factor = "GFBio", .before = 1)
  ),
  "Summaries/chisq_posthoc.txt",
  quote = F,
  sep = "\t",
  row.names = F
)

write.table(
  bind_rows(
    chisq.entry.db.oq.posthoc %>% mutate(factor = "entry DB", .before = 1),
    chisq.mixs.oq.posthoc %>% mutate(factor = "MIxS", .before = 1),
    chisq.gfbio.oq.posthoc %>% mutate(factor = "GFBio", .before = 1)
  ),
  "Summaries/chisq_overall_quality.txt",
  quote = F,
  sep = "\t",
  row.names = F
)


### how many studies / runs submitted with primers ####

require(tidyverse)
require(reshape)

# case study Riedinger et al. 2024
# V3V4 16S amplicons
# manually confirmed primer presence
primer_present <- read.table("Case_studies/primer_presense_start_end_DR.csv", h = T, sep = ",") %>% 
  filter(primer_start_end != "other/unclear")

# determined if paired ends already merged prior to submission
# taken from: https://git.io-warnemuende.de/riedinge/datamining_vulnificus_DR_CH/src/branch/main/snakemake_wf/assets
#   pe.separate <- scan(what = "character", file = "Case_studies/sample_list_2.txt", sep = "\n")
#   pe.merged <- scan(what = "character", file = "Case_studies/sample_list_1.txt", sep = "\n")
# there was a logical error in the assumption about data processing at ENA in the above selection
# Any run with 3 files had uploaded 2 files to begin with --> the third file is generated by ENA for unmatched pairs
# Runs with 1 file are either merged or all reads were unpaired in submission with 2 files (e.g. wrong files uploaded because of mix-up in file name)

# summarize metadata for those runs
ena.assess.primer <- ena.assess[match(primer_present$run_accession, ena.assess$run_accession), ] %>% 
  mutate(
    primer_present = as.factor(ifelse(primer_present$primer_start_end == "yes V3_V4", "included", "removed")),
    merged_status = as.factor(ifelse(sapply(strsplit(fastq_ftp, ";"), length) == 1, "pe_merged", "pe_separate"))
  )
all.equal(primer_present$run_accession, ena.assess.primer$run_accession)
primer.study <- ena.assess.primer %>% select(study_accession, mixs, entry_db, gfbio, primer_present, merged_status) %>% unique()
primer.study %>% pull(study_accession) %>% table() %>% table()
# submission strategy consistent per study
# perform further assessment at study level

stats.primer.usecase <- primer.study %>% 
  group_by(entry_db, mixs, gfbio, primer_present, merged_status, .drop = F) %>% 
  summarize(n_data = n())
write.table(stats.primer.usecase, "Summaries/summary_primer_usecase.txt", sep = "\t", quote = F, row.names = F)
saveRDS(stats.primer.usecase, "Summaries/summary_primer_usecase.rds")
write.table(primer.study, "Summaries/riedinger_usecase.txt", sep = "\t", quote = F, row.names = F)
saveRDS(primer.study, "Summaries/riedinger_usecase.rds")

cast(stats.primer.usecase, "entry_db + primer_present + merged_status ~ .", value = "n_data", fun.aggregate = sum) %>% 
  clean_names() %>% 
  group_by(entry_db) %>% 
  mutate(percent = round(all/sum(all) * 100, 2))
# only NCBI with already merged reads
cast(stats.primer.usecase, "mixs + primer_present + merged_status ~ .", value = "n_data", fun.aggregate = sum) %>% 
  clean_names() %>% 
  group_by(mixs) %>% 
  mutate(percent = round(all/sum(all) * 100, 2))
# higher percentage of primer removed unmerged
cast(stats.primer.usecase, "gfbio + primer_present + merged_status ~ .", value = "n_data", fun.aggregate = sum) %>% 
  clean_names() %>% 
  group_by(gfbio) %>% 
  mutate(percent = round(all/sum(all) * 100, 2))
# too few studies with gfbio to draw conclusions