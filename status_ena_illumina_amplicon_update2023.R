# Assessing status of metadata on ENA for Illumina amplicon studies 
# from metagenomic material (taxid: ecological metagenomes)

### prepare environment and read data ####

setwd("C:/Users/chassenrueck/Documents/Additional_projects/Data_mining_commentary/Update_2023/")
require(data.table)
require(reshape)
require(dplyr)
require(tidyverse)
require(janitor)
require(furrr)
require(dlfUtils)
require(R.utils)

source("C:/Users/chassenrueck/Documents/Repos/ENA-metadata/plot_ena_summaries.R")

# save.image("status_ena_illumina_amplicon_update2023.Rdata")
load("status_ena_illumina_amplicon_update2023.Rdata")

# read ENA search results
ena.out <- fread(
  "results_read_run_tsv.txt",
  h = T,
  sep = "\t",
  quote = ""
)

# which columns are completely empty
empty.metadata <- colnames(ena.out)[apply(ena.out, 2, function(x) sum(is.na(x) | x == "") == length(x))]
length(empty.metadata)

# remove empty columns
ena.out <- ena.out[, (empty.metadata):=NULL]
dim(ena.out)
gc()
# 891315 entries, 113 fields with data

# studies with publication
xref.out <- read.table("xref_out/studies_with_paper_nr.txt", h = F, sep = "\t")

# INSDC country list
insdc_countries <- scan("insdc_country_list.txt", sep = "\n", what = "character")

# parameters of interest for analysis
select_params <- c(
  "run_accession",
  "sample_accession",
  "altitude",
  "base_count",
  "broad_scale_environmental_context",
  "broker_name",
  "center_name",
  "checklist",
  "collection_date",
  "collection_date_end",
  "collection_date_start",         
  "country",
  "depth",
  "description",
  "elevation",
  "environment_biome",
  "environment_feature",
  "environment_material",
  "environmental_medium",
  "environmental_sample",
  "experiment_accession",
  "experiment_alias",
  "experiment_title",
  "fastq_bytes",
  "first_created",
  "first_public",
  "instrument_model",
  "isolation_source",
  "last_updated",
  "lat",
  "library_construction_protocol",
  "library_name",
  "library_selection",
  "local_environmental_context",
  "location",
  "location_end",
  "location_start",
  "lon",
  "ncbi_reporting_standard",
  "nominal_length",
  "study_accession",
  "project_name",
  "read_count",
  "run_alias",
  "sample_alias",
  "sample_description",
  "sample_title",
  "tax_id",
  "scientific_name",                 
  "secondary_sample_accession",
  "secondary_study_accession",
  "study_alias",
  "study_title",
  "target_gene"
)
ena.out <- ena.out[, ..select_params]
dim(ena.out)
gc()
# 891315 entries, 54 data fields

# define evaluation criteria (grouping categories)
mixs.checklists <- paste0("ERC0000", c(12:25, 31))
ena.assess <- data.frame(
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  MIxS = grepl("MI", ena.out$ncbi_reporting_standard, fixed = T) | ena.out$checklist %in% mixs.checklists,
  study.accnos = ena.out$study_accession,
  sample.accnos = ena.out$sample_accession,
  exp.accnos = ena.out$experiment_accession
)
ena.assess$year.created <- as.factor(ena.assess$year.created)
apply(ena.assess[, c("year.created", "gfbio", "MIxS")], 2, table)


### parse data from xml ####

# run this on bio-48 with 32 parallel threads in furrr
# conda activate r-4.2.2
# cd /home/hassenru/ENA_data_mining/xml_parsed_csv
plan(multicore, workers = 32)
location_xml <- scan("sample_attributes_location_xy.txt", what = "character", sep = "\n")
location_z_xml <- scan("sample_attributes_location_z.txt", what = "character", sep = "\n")
location_region_xml <- scan("sample_attributes_location_region.txt", what = "character", sep = "\n")
primers_xml <- scan("sample_attributes_primers.txt", what = "character", sep = "\n")
date_xml <- scan("sample_attributes_date.txt", what = "character", sep = "\n")
xml_list <- list.files(path = "sample_csv", pattern = "*.csv")

xml_df <- future_map(
  xml_list,
  function(x) {
    tmp <- read.table(
      paste0("sample_csv/", x),
      h = T,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = F,
      colClasses = "character"
    )
  }
)

df_loc <- future_map_dfr(
  xml_df,
  function(x) {
    x[, colnames(x) %in% c("PRIMARY ID", location_xml), drop = F]
  }
)
write.table(df_loc, "df_parsed_location.txt", sep = "\t", quote = F, row.names = F)
df_z <- future_map_dfr(
  xml_df,
  function(x) {
    x[, colnames(x) %in% c("PRIMARY ID", location_z_xml), drop = F]
  }
)
write.table(df_z, "df_parsed_z.txt", sep = "\t", quote = F, row.names = F)
df_region <- future_map_dfr(
  xml_df,
  function(x) {
    x[, colnames(x) %in% c("PRIMARY ID", location_region_xml), drop = F]
  }
)
write.table(df_region, "df_parsed_region.txt", sep = "\t", quote = F, row.names = F)
df_primers <- future_map_dfr(
  xml_df,
  function(x) {
    x[, colnames(x) %in% c("PRIMARY ID", primers_xml), drop = F]
  }
)
write.table(df_primers, "df_parsed_primers.txt", sep = "\t", quote = F, row.names = F)
df_date <- future_map_dfr(
  xml_df,
  function(x) {
    x[, colnames(x) %in% c("PRIMARY ID", date_xml), drop = F]
  }
)
write.table(df_date, "df_parsed_date.txt", sep = "\t", quote = F, row.names = F)
# save.image("collect_xml_info.Rdata")

# move to local computer
df_xml_xy <- data.frame(
  fread("xml_parsed_csv/df_parsed_location.txt", h = T,  sep = "\t",  quote = ""),
  check.names = F,
  row.names = 1
)
# df_xml_z <- data.frame(
#   fread("xml_parsed_csv/df_parsed_z.txt", h = T,  sep = "\t",  quote = ""),
#   check.names = F,
#   row.names = 1
# )
df_xml_region <- data.frame(
  fread("xml_parsed_csv/df_parsed_region.txt", h = T,  sep = "\t",  quote = ""),
  check.names = F,
  row.names = 1
)
df_xml_primers <- data.frame(
  fread("xml_parsed_csv/df_parsed_primers.txt", h = T,  sep = "\t",  quote = ""),
  check.names = F,
  row.names = 1
)
df_xml_date <- data.frame(
  fread("xml_parsed_csv/df_parsed_date.txt", h = T,  sep = "\t",  quote = ""),
  check.names = F,
  row.names = 1
)


# check xy data (longitude, latitude) ####
# remove samples, which already have values in TSV output
df_xml_xy <- df_xml_xy[rownames(df_xml_xy) %in% ena.out$sample_accession[is.na(ena.out$lat) | is.na(ena.out$lon)], ]
# remove all NA entries from parsed xml table
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# inspect and clean individual parameters
unique(sort(df_xml_xy$`lat lon`))
df_xml_xy$`lat lon`[
  df_xml_xy$`lat lon` %in% c("-", "~50S", "0", "not available", "Not available", "not collected", "Not collected", "Not Collected", "not pplicable", "Not Recorded", "restricted access", "restricted_access")
] <- NA
tmp_accnos <- rownames(df_xml_xy)[!is.na(df_xml_xy$`lat lon`)]
# remove sample with validated information
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$lat_lon))
df_xml_xy$lat_lon[
  df_xml_xy$lat_lon %in% c("missing", "Missing", "N/A", "na", "not applicable", "Not applicable", "Not Applicable", "not collected", "Not collected", "not provided", "restricted access")
] <- NA
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$lat_lon)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter (assuming that those 2 go together)
unique(sort(df_xml_xy$`geographic location (longitude)`))
unique(sort(df_xml_xy$`geographic location (latitude)`))
df_xml_xy$`geographic location (longitude)`[
  df_xml_xy$`geographic location (longitude)` %in% c("not collected", "not provided", "restricted access", "1a") |
    df_xml_xy$`geographic location (latitude)` %in% c("not collected", "not provided", "restricted access", "1a")
] <- NA
df_xml_xy$`geographic location (latitude)`[
  df_xml_xy$`geographic location (longitude)` %in% c("not collected", "not provided", "restricted access", "1a") |
    df_xml_xy$`geographic location (latitude)` %in% c("not collected", "not provided", "restricted access", "1a")
] <- NA
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$`geographic location (longitude)`) & !is.na(df_xml_xy$`geographic location (latitude)`)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
# geographic location: country information
tmp <- full_join(
  df_xml_region %>% rownames_to_column(),
  df_xml_xy %>% rownames_to_column() %>% select(rowname, 'geographic location'),
  by = "rowname"
)
df_xml_region <- tmp
rm(tmp)
df_xml_xy$`geographic location` <- NULL
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$`latitude and longitude`))
df_xml_xy$`latitude and longitude` <- NULL
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter (sampling event all NA)
df_xml_xy$`sampling event, latitude, start` <- NULL
df_xml_xy$`sampling event, latitude, end` <- NULL
df_xml_xy$`sampling event, longitude, start` <- NULL
df_xml_xy$`sampling event, longitude, end` <- NULL
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$`geographic location (latitudeandlongitude)`))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$`geographic location (latitudeandlongitude)`)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$`lat_lon Site A`))
unique(sort(df_xml_xy$`lat_lon Site B`))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$`lat_lon Site A`) & !is.na(df_xml_xy$`lat_lon Site B`)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$mixed_lat_lon))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$mixed_lat_lon)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$`lat_lo N`))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$'lat_lo N')])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$`Geographic Location`))
# Geographic Location: country information
tmp <- full_join(
  df_xml_region,
  df_xml_xy %>% rownames_to_column() %>% select(rowname, 'Geographic Location'),
  by = "rowname"
)
df_xml_region <- tmp
rm(tmp)
df_xml_xy$`Geographic Location` <- NULL
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$latitude_and_longitude))
df_xml_xy$latitude_and_longitude[
  df_xml_xy$latitude_and_longitude %in% c("not applicable", "not collected")
] <- NA
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$latitude_and_longitude)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameters
df_xml_xy$gps_number <- NULL
df_xml_xy$geographic_location <- NULL
df_xml_xy$Geo <- NULL
df_xml_xy$`spatial position` <- NULL
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameter
unique(sort(df_xml_xy$lat_lon_range))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$lat_lon_range)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
# next parameters
unique(sort(df_xml_xy$Latitude))
unique(sort(df_xml_xy$Longitude))
tmp_accnos <- c(tmp_accnos, rownames(df_xml_xy)[!is.na(df_xml_xy$Latitude) & !is.na(df_xml_xy$Longitude)])
df_xml_xy <- df_xml_xy[!rownames(df_xml_xy) %in% tmp_accnos, ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[!apply(df_xml_xy, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_xy <- df_xml_xy[, !apply(df_xml_xy, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
rm(df_xml_xy)

# check region data (for any additional xy info) ####
df_xml_region <- df_xml_region %>% column_to_rownames("rowname")
df_xml_region <- df_xml_region[rownames(df_xml_region) %in% ena.out$sample_accession[ena.out$country == ""], ]
df_xml_region <- df_xml_region[!apply(df_xml_region, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_region <- df_xml_region[, !apply(df_xml_region, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
unique(sort(df_xml_region$`geo loc name`))
df_xml_region$`geo loc name`[
  df_xml_region$`geo loc name` %in% c("not available", "Not available", "not collected", "restricted access")
] <- NA
unique(sort(df_xml_region$geo_loc_name))
df_xml_region$geo_loc_name[
  df_xml_region$geo_loc_name %in% c("missing", "not applicable", "Not applicable", "not collected", "not provided")
] <- NA
df_xml_region$`geographic location (region and locality)` <- NULL
unique(sort(df_xml_region$`geographic location (country and/or sea)`))
df_xml_region$`geographic location (country and/or sea)`[
  df_xml_region$`geographic location (country and/or sea)` %in% c("not collected")
] <- NA
df_xml_region$country <- NULL
unique(sort(df_xml_region$`geographic location (countryand/orsea,region)`))
df_xml_region$`sampling event, device collection area` <- NULL
unique(sort(df_xml_region$region))
tmp <- full_join(
  df_xml_primers %>% rownames_to_column(),
  df_xml_region %>% rownames_to_column() %>% select(rowname, region),
  by = "rowname"
)
df_xml_primers <- tmp %>% column_to_rownames("rowname")
rm(tmp)
df_xml_region$region <- NULL
df_xml_region$location <- NULL
df_xml_region$`sample location` <- NULL
df_xml_region$continent <- NULL
unique(sort(df_xml_region$geolocname))
unique(sort(df_xml_region$`collection site`))
df_xml_region$`location code` <- NULL
df_xml_region$geo_region <- NULL
df_xml_region$`Location ID` <- NULL
df_xml_region$`sample_name,sample_title,bioproject_accession,organism,host,isolation_source,collection_date,geo_loc_name,lat_lon` <- NULL
df_xml_region$city_ID <- NULL
df_xml_region$Collection_area <- NULL
df_xml_region$SampleOrigin <- NULL
df_xml_region$sample_origin <- NULL
df_xml_region$`collection point` <- NULL
df_xml_region$Location <- NULL
df_xml_region$collection_site_id <- NULL
df_xml_region$country_of_origin <- NULL
df_xml_region$locality_name <- NULL
df_xml_region$coll_site_geo_feat <- NULL
df_xml_region$Environtal_location <- NULL
df_xml_region$`geographic location` <- NULL
df_xml_region$`Geographic Location` <- NULL
# no hidden xy information
accnos_xml_xy <- tmp_accnos
df_xml_region <- df_xml_region[!apply(df_xml_region, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_region <- df_xml_region[, !apply(df_xml_region, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
accnos_xml_region <- rownames(df_xml_region)

# check primer info ####
df_xml_primers <- df_xml_primers[rownames(df_xml_date) %in% ena.out$sample_accession[ena.out$target_gene == ""], ]
df_xml_primers <- df_xml_primers[!apply(df_xml_primers, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_primers <- df_xml_primers[, !apply(df_xml_primers, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
unique(sort(df_xml_primers$`target gene`))
unique(sort(df_xml_primers$`pcr primers`))
unique(sort(df_xml_primers$`target subfragment`))
unique(sort(df_xml_primers$`target gene.1`))
df_xml_primers$seq_lib <- NULL
unique(sort(df_xml_primers$primer_r))
unique(sort(df_xml_primers$primer_r_reference))
unique(sort(df_xml_primers$primer_f))
unique(sort(df_xml_primers$primer_f_reference))
df_xml_primers$barcodesequence <- NULL
df_xml_primers$PrimerNumber <- NULL
df_xml_primers$DNA <- NULL
unique(sort(df_xml_primers$domain))
unique(sort(df_xml_primers$locus_tag))
df_xml_primers$`nucleic acid amplification` <- NULL
unique(sort(df_xml_primers$`Forward_MT-FS`))
df_xml_primers$Sample_Barcode_Name <- NULL
unique(sort(df_xml_primers$`Reverse_MT-FS`))
df_xml_primers$Sample_Barcode_Sequence <- NULL
unique(sort(df_xml_primers$`pcr primers.1`))
df_xml_primers$`amplicon pool` <- NULL
df_xml_primers$`barcode sequence` <- NULL
unique(sort(df_xml_primers$rrna_region))
unique(sort(df_xml_primers$targetgene))
table(nchar(df_xml_primers$`reverse primer`))
df_xml_primers$`reverse primer`[
  !df_xml_primers$`reverse primer` %in% c("bact_805R") & nchar(df_xml_primers$`reverse primer`) < 18
] <- NA
unique(sort(df_xml_primers$`forward primer`))
df_xml_primers$`forward primer`[
  !df_xml_primers$`forward primer` %in% c("bact_341F") & nchar(df_xml_primers$`forward primer`) < 17
] <- NA
df_xml_primers$barcode <- NULL
unique(sort(df_xml_primers$`marker gene`))
df_xml_primers$its.library <- NULL
df_xml_primers$ssu_dna <- NULL
df_xml_primers$its_dna <- NULL
df_xml_primers$its_g <- NULL
df_xml_primers$ssu_g <- NULL
df_xml_primers$ssu.library <- NULL
df_xml_primers$`library construction method` <- NULL
unique(sort(df_xml_primers$forwardprimersequence))
unique(sort(df_xml_primers$barcodenames))
unique(sort(df_xml_primers$reverseprimersequence))
unique(sort(df_xml_primers$primer))
df_xml_primers$primer[
  grepl("EMP", df_xml_primers$primer) | df_xml_primers$primer %in% c("f", "cons", "forward", "reverse", "Rust")
] <- NA
unique(sort(df_xml_primers$library))
df_xml_primers$library[
  !grepl("_V4_|_V34_", df_xml_primers$library)
] <- NA
df_xml_primers$reversebarcodesequence <- NULL
df_xml_primers$forwardbarcodesequence <- NULL
unique(sort(df_xml_primers$amplicon))
df_xml_primers$`genetic material` <- NULL
unique(sort(df_xml_primers$primers))
df_xml_primers$probe <- NULL
unique(sort(df_xml_primers$`linker primer sequence`))
unique(sort(df_xml_primers$`barcode name`))
unique(sort(df_xml_primers$rRNA_fragment))
df_xml_primers$Forward_Reverse_Sequence <- NULL
unique(sort(df_xml_primers$reverseprimer))
unique(sort(df_xml_primers$gene))
df_xml_primers$`Concatenated Barcode` <- NULL
df_xml_primers$`forward barcode` <- NULL
df_xml_primers$`barcode id` <- NULL
df_xml_primers$`Barcode + Linker UniTag2` <- NULL
df_xml_primers$`Barcode + Linker UniTag1` <- NULL
df_xml_primers$`Linker UniTag1` <- NULL
df_xml_primers$`Linker UniTag2` <- NULL
unique(sort(df_xml_primers$`molecule type`))
df_xml_primers$template <- NULL
unique(sort(df_xml_primers$`Reverse primer sequence (5' -> '3)`))
unique(sort(df_xml_primers$`Forwar primer ID`))
df_xml_primers$`DNA sample type` <- NULL
unique(sort(df_xml_primers$`Forward primer sequence (5' -> '3)`))
unique(sort(df_xml_primers$`Reverse primer ID`))
unique(sort(df_xml_primers$target_region))
df_xml_primers$`barcode number` <- NULL
df_xml_primers$`For-rev` <- NULL
unique(sort(df_xml_primers$target))
df_xml_primers$target[
  df_xml_primers$target %in% c("DNA")
] <- NA
df_xml_primers$`Genetic material` <- NULL
unique(sort(df_xml_primers$`gene name`))
unique(sort(df_xml_primers$`target community`))
df_xml_primers$rev_barcode <- NULL
df_xml_primers$fwd_barcode <- NULL
df_xml_primers$fwd_rev_read <- NULL
unique(sort(df_xml_primers$`primers for nifH amplification`))
df_xml_primers$amplicon_library_type <- NULL
unique(sort(df_xml_primers$amplicon_sequenced))
unique(sort(df_xml_primers$loci))
unique(sort(df_xml_primers$`rev primer`))
unique(sort(df_xml_primers$locus))
df_xml_primers$locus[
  df_xml_primers$locus %in% as.character(1:9)
] <- NA
unique(sort(df_xml_primers$Fwd_primer))
unique(sort(df_xml_primers$`16s region`))
unique(sort(df_xml_primers$target_domain))
unique(sort(df_xml_primers$primers_target))
unique(sort(df_xml_primers$`primer set`))
unique(sort(df_xml_primers$Fwd_primer_seq))
unique(sort(df_xml_primers$Rev_primer_name))
unique(sort(df_xml_primers$Rev_primer_seq))
unique(sort(df_xml_primers$Fwd_primer_name))
unique(sort(df_xml_primers$`primer pair`))
df_xml_primers$dnatag <- NULL
unique(sort(df_xml_primers$`sample_name,organism,host,isolation_source,collection_date,geo_loc_name,lat_lon`))
unique(sort(df_xml_primers$`marker used`))
df_xml_primers$`Collection Date Primer 1 (EUK14)` <- NULL
df_xml_primers$`Collection Date Primers 2 (EUK15 + DIV4)` <- NULL
df_xml_primers$`Replicate Primers 2` <- NULL
df_xml_primers$`Replicate Primer 1` <- NULL
df_xml_primers$primerID <- NULL
df_xml_primers$primer_combination <- NULL
df_xml_primers$barcode_combination <- NULL
unique(sort(df_xml_primers$sequence))
df_xml_primers$sequence[
  df_xml_primers$sequence %in% c("forward")
] <- NA
unique(sort(df_xml_primers$Locus_target))
unique(sort(df_xml_primers$MolecularMarker))
df_xml_primers$molecule <- NULL
unique(sort(df_xml_primers$gene_region))
unique(sort(df_xml_primers$`sample_name,sample_title,bioproject_accession,organism,host,isolation_source,collection_date,geo_loc_name,lat_lon`))
unique(sort(df_xml_primers$`library construction protocol`))
df_xml_primers$linker <- NULL
df_xml_primers$primer_id <- NULL
unique(sort(df_xml_primers$rRNA))
unique(sort(df_xml_primers$`gene target`))
unique(sort(df_xml_primers$Targeted))
unique(sort(df_xml_primers$`sequence type`))
df_xml_primers$Fwd_Primer_Barcode <- NULL
df_xml_primers$Barcodes <- NULL # not really sure... could maybe be primers...
df_xml_primers$Rev_Primer_Barcode <- NULL
df_xml_primers$PCR_product <- NULL
unique(sort(df_xml_primers$`Reverse Primer (926R)`))
unique(sort(df_xml_primers$`Forward Primer (515F)`))
unique(sort(df_xml_primers$`Forward Primer (H279)`))
# I am not checking the other parameters that follow a similar nomenclature
df_xml_primers$template_for_pcr <- NULL
unique(sort(df_xml_primers$`Forward Primer (TAReuk454FWD1)`))
unique(sort(df_xml_primers$Amplicon_marker))
unique(sort(df_xml_primers$amplified_dna_region))
unique(sort(df_xml_primers$Gene_hypervariable_regions))
unique(sort(df_xml_primers$Gene_sequenced))
df_xml_primers$Gene_sequenced[
  df_xml_primers$Gene_sequenced %in% c("Metagenome")
] <- NA
unique(sort(df_xml_primers$primerset))
df_xml_primers$`DNA_sample.id;unique_sampID;treatment;myc_assoc;incubated;sample_type;type;isolate;mel_status;incub_period.days` <- NULL
df_xml_primers$amplicon_sequencing_MBL_id <- NULL
unique(sort(df_xml_primers$sequenced_region))
df_xml_primers$linkerprimer <- NULL
df_xml_primers$primer_number <- NULL
unique(sort(df_xml_primers$primer_pairs))
unique(sort(df_xml_primers$fwd_primer))
unique(sort(df_xml_primers$Primerset))
df_xml_primers$`samplecode(markerregionvineyardmgmtyearseason)` <- NULL
df_xml_primers$`SampleCode (Marker_Region_Vineyard_Mgmt_Year_Season)` <- NULL
df_xml_primers$`SampleCode (Region_Vineyard_Mgmt_Year_Season)` <- NULL
unique(sort(df_xml_primers$`reverse primer sequence`))
unique(sort(df_xml_primers$fwd.primer.seq))
unique(sort(df_xml_primers$rev.primer.seq))
unique(sort(df_xml_primers$primer_set_ITS))
unique(sort(df_xml_primers$primer_set_28S))
df_xml_primers$nucl_acid_amp <- NULL
df_xml_primers$'16SrRNAgene_libraries_linked_to_sample' <- NULL
df_xml_primers$barcode_library_1 <- NULL
df_xml_primers$sequences <- NULL
df_xml_primers$'sample barcode' <- NULL
df_xml_primers$'Reverse_Index_ Read_I1' <- NULL
df_xml_primers$'Forward_Index_Read_I2' <- NULL
df_xml_primers$PrimerPosition <- NULL
df_xml_primers$`eDNA treatment` <- NULL
df_xml_primers$rRNA_16S <- NULL
# colnames(df_xml_primers)[114]
# unique(sort(df_xml_primers[, 114]))
df_xml_primers$Target_PCT_Name <- NULL
df_xml_primers$Target_PCT_No <- NULL
unique(sort(df_xml_primers$`Marker type`))
unique(sort(df_xml_primers$marker_region))
unique(sort(df_xml_primers$`Amplified Region`))
unique(sort(df_xml_primers$LinkerPrimerSequence))
df_xml_primers$BarcodeSequence <- NULL
df_xml_primers$BarcodeName <- NULL
unique(sort(df_xml_primers$pcr_primers))
unique(sort(df_xml_primers$target_gene))
unique(sort(df_xml_primers$target_subfragment))
unique(sort(df_xml_primers$Amplicon_type))
unique(sort(df_xml_primers$Barcode))
unique(sort(df_xml_primers$forward_primer))
unique(sort(df_xml_primers$reverse_primer))
unique(sort(df_xml_primers$Domain))
unique(sort(df_xml_primers$ww_surv_target_1_gene))
df_xml_primers$ww_surv_target_1 <- NULL
df_xml_primers$ww_surv_target_2 <- NULL
unique(sort(df_xml_primers$`ITS/16S`))
unique(sort(df_xml_primers$amplicon_type))
unique(sort(df_xml_primers$Primer))
unique(sort(df_xml_primers$PrimerDesc))
unique(sort(df_xml_primers$ReversePrimer))
unique(sort(df_xml_primers$`16S/ITS`)) # not all parameters using this nomenclature checked
df_xml_primers$BARCODE <- NULL
unique(sort(df_xml_primers$Target_gene))
unique(sort(df_xml_primers$`Pcr_primer pair`))
unique(sort(df_xml_primers$Target_subfragment))
unique(sort(df_xml_primers$'16s_region'))
df_xml_primers$`barcode sequences` <- NULL
unique(sort(df_xml_primers$AmpliconType))
df_xml_primers$Library <- NULL
unique(sort(df_xml_primers$Marker))
unique(sort(df_xml_primers$ForwardPrimerSequence))
df_xml_primers$barcode_forward <- NULL
df_xml_primers$barcode_reverse <- NULL
unique(sort(df_xml_primers$`Primer set`))
unique(sort(df_xml_primers$Gene))
df_xml_primers$DNA_type_sequenced <- NULL
unique(sort(df_xml_primers$`Reverse Primer`))
unique(sort(df_xml_primers$`Forward Primer`))
unique(sort(df_xml_primers$Sequencing_Marker))
unique(sort(df_xml_primers$`Target gene`))
unique(sort(df_xml_primers$`dna region`))
unique(sort(df_xml_primers$Primers))
df_xml_primers$Sequence <- NULL
df_xml_primers$edna_label <- NULL
df_xml_primers$`#template` <- NULL
df_xml_primers <- df_xml_primers[!apply(df_xml_primers, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_primers <- df_xml_primers[, !apply(df_xml_primers, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
accnos_xml_primers <- rownames(df_xml_primers)


# check collection date ####
df_xml_date <- df_xml_date[rownames(df_xml_date) %in% ena.out$sample_accession[(is.na(ena.out$collection_date) | ena.out$collection_date == "") & (ena.out$collection_date_start == "" | ena.out$collection_date_end == "")], ]
df_xml_date <- df_xml_date[!apply(df_xml_date, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_date <- df_xml_date[, !apply(df_xml_date, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
unique(sort(df_xml_date$`collection date`))
df_xml_date$`collection date`[
  df_xml_date$`collection date` %in% c("missing: control sample", "not collected")
] <- NA
unique(sort(df_xml_date$collection_date))
df_xml_date$collection_date[
  df_xml_date$collection_date %in% c("missing", "not applicable", "Not applicable", "not collected")
] <- NA
unique(sort(df_xml_date$`event date/time start`))
df_xml_date$time <- NULL
df_xml_date$timepoint <- NULL
unique(sort(df_xml_date$`Event Date/Time Start`))
df_xml_date$`sampling event, date/time, end` <- NULL
df_xml_date$`sampling event, date-time, start` <- NULL
unique(sort(df_xml_date$date))
unique(sort(df_xml_date$year))
df_xml_date$day_sample <- NULL
df_xml_date$`sampling day` <- NULL
unique(sort(df_xml_date$collection_timestamp))
df_xml_date$`time point` <- NULL
unique(sort(df_xml_date$collectiondate))
unique(sort(df_xml_date$mixed_collection_date))
df_xml_date$`Collection Date` <- NULL
unique(sort(df_xml_date$Collection_date))
df_xml_date$Collection_date[
  df_xml_date$Collection_date %in% c("na")
] <- NA
df_xml_date$incubation_date <- NULL
df_xml_date$Months <- NULL
df_xml_date$sample_collection_week <- NULL
unique(sort(df_xml_date$location_survey_date))
df_xml_date$extraction_date <- NULL
df_xml_date$Time_days <- NULL
df_xml_date$Extraction_Date <- NULL
df_xml_date$DNA_Allocation_Date <- NULL
df_xml_date$site_date <- NULL
df_xml_date$Time <- NULL
df_xml_date$`extraction date` <- NULL
df_xml_date <- df_xml_date[!apply(df_xml_date, 1, function(x) sum(is.na(x)) == length(x)), ] # repeat after each cleaning step
df_xml_date <- df_xml_date[, !apply(df_xml_date, 2, function(x) sum(is.na(x)) == length(x))] # repeat after each cleaning step
accnos_xml_date <- rownames(df_xml_date)


# remove xml output to reduce workspace size ####
rm(df_xml_date, df_xml_primers, df_xml_region, df_xml_z)


### basic stats ####

nrow(ena.assess)
# 891315
length(unique(ena.assess$study.accnos))
# 15214
length(unique(ena.assess$sample.accnos))
# 752905
length(unique(ena.assess$exp.accnos))
# 888976
length(unique(ena.assess$sample.accnos)) - 746689 # number of parsable XML with data
# 6216
summary(c(table(ena.assess$sample.accnos)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.184   1.000 776.000 

# paired end data should have 2 files per run
table(sapply(strsplit(ena.out$fastq_bytes, ";"), function(x) length(x)))
# 0      1      2      3 
# 14152  77960 785157  14046 
# 12% of data sets not archived correctly

sum(unlist(sapply(strsplit(ena.out$fastq_bytes, ";"), function(x) as.numeric(x))))



### check taxid ####
sum(is.na(ena.out$tax_id))
sum(ena.out$tax_id == "")
# taxid is always provided


### check base count ####
sum(is.na(ena.out$base_count))
sum(ena.out$base_count == "")
sum(ena.out$base_count == 0)
# 14145 runs without data
ena.assess$no.data <- ifelse(ena.out$base_count == 0, "no_data", "with_data")
ena.assess$no.data <- factor(ena.assess$no.data, levels = c("with_data", "no_data"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + no.data ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_no_data.txt", sep = "\t", quote = F, row.names = F)
ena.assess %>% filter(gfbio & MIxS & no.data == "no_data") %>% View()
# this could be an internal file processing error or orphan projects


### check collection date ####
sum(is.na(ena.out$collection_date))
sum(ena.out$collection_date == "", na.rm = T)
sum(is.na(ena.out$collection_date_end))
sum(ena.out$collection_date_end == "", na.rm = T)
sum(is.na(ena.out$collection_date_start))
sum(ena.out$collection_date_start == "", na.rm = T)
# check syntax
sum(
  grepl(
    "(^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)",
    ena.out$collection_date
  )
)
# only concrete numbers
sum(
  grepl(
    "(^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$)",
    ena.out$collection_date
  )
)

ena.assess$collection.date.available <- ifelse(
  !(is.na(ena.out$collection_date) | ena.out$collection_date == ""),
  "yes",
  ifelse(
    !(ena.out$collection_date_start == "" | ena.out$collection_date_end == ""),
    "period",
    ifelse(
      ena.out$sample_accession %in% accnos_xml_date,
      "xml",
      "no"
    )
  )
)
ena.assess$collection.date.available <- factor(ena.assess$collection.date.available, levels = c("yes", "period", "xml", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + collection.date.available ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_collection_date_available.txt", sep = "\t", quote = F, row.names = F)
ena.assess %>% filter(gfbio & !MIxS & collection.date.available == "no") %>% View()
# could also be because of negative controls or experimental settings


### check lat/lon ####
sum(is.na(ena.out$lat))
sum(is.na(ena.out$lon))
sum(is.na(ena.out$lon) & !is.na(ena.out$lat))
# 202024 NAs
sum(ena.out$lat == "", na.rm = T)
sum(ena.out$lon == "", na.rm = T)
# no empty cells
sum(is.na(ena.out$location))
sum(ena.out$location == "", na.rm = T)
sum(is.na(ena.out$location_start))
sum(ena.out$location_start == "", na.rm = T)
sum(is.na(ena.out$location_end))
sum(ena.out$location_end == "", na.rm = T)
sum(is.na(ena.out$lon) & ena.out$location != "")
# there are no cases where location* has values, but lat/lon don't

# check syntax
sum(
  grepl(
    "(^[+-]?[0-9]+.?[0-9]{0,8}$)|(^not collected$)|(^not provided$)|(^restricted access$)|(^missing: control sample$)|(^missing: sample group$)|(^missing: synthetic construct$)|(^missing: lab stock$)|(^missing: third party data$)|(^missing: data agreement established pre-2023$)|(^missing: endangered species$)|(^missing: human-identifiable$)",
    ena.out$lon
  )
)
# only concrete numbers
sum(
  grepl(
    "(^[+-]?[0-9]+.?[0-9]{0,8}$)",
    ena.out$lon
  )
)
# no cases where non-numbers are inlcuded in ena TSV output

ena.assess$lat_lon <- ifelse(
  !is.na(ena.out$lat) & !is.na(ena.out$lon),
  "yes", 
  ifelse(
    ena.out$sample_accession %in% accnos_xml_xy,
    "xml",
    "no"
  )
)
ena.assess$lat_lon <- factor(ena.assess$lat_lon, levels = c("yes", "xml", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + lat_lon ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_lat_lon.txt", sep = "\t", quote = F, row.names = F)


### check country information and syntax ####
sum(is.na(ena.out$country))
# no NAs
sum(ena.out$country == "", na.rm = T)
# 107242 empty cells
# search for country/sea name anywhere in entry
test <- sapply(insdc_countries, function(x) grepl(x, ena.out$country, fixed = T))
ena.assess$country.available <- ifelse(
  ena.out$sample_accession %in% accnos_xml_region,
  "xml",
  ifelse(
    ena.out$country == "",
    "no",
    ifelse(
      apply(test, 1, any),
      "insdc_country",
      "other"
    )
  )
)
ena.assess$country.available <- factor(ena.assess$country.available, levels = c("insdc_country", "other", "xml", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + country.available ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_country_available.txt", sep = "\t", quote = F, row.names = F)
# ENA also allows: missing (for various reasons, e.g. control samples), not applicable, not collected, not provided, restricted access
# compared to previous version of API, controlled vocabulary is no used


### check z-coordinates ####
sum(is.na(ena.out$altitude))
sum(is.na(ena.out$depth))
sum(is.na(ena.out$elevation))
# 881119, 770, 853387 NAs
sum(ena.out$altitude == "", na.rm = T) # no empty cells
sum(ena.out$depth == "", na.rm = T) # 659294 empty cells
sum(ena.out$elevation == "", na.rm = T) # no empty cells
# hardly ever provided --> at least 1 value should be available
# regular expression required
# altitude:        (0|((0\.)|([1-9][0-9]*\.?))[0-9]*)([Ee][+-]?[0-9]+)? 
# depth:           (0|((0\.)|([1-9][0-9]*\.?))[0-9]*)([Ee][+-]?[0-9]+)? 
# elevation": [+-]?(0|((0\.)|([1-9][0-9]*\.?))[0-9]*)([Ee][+-]?[0-9]+)? 
test_pa <- data.frame(
  altitude = !(is.na(ena.out$altitude) | ena.out$altitude == ""),
  depth = !(is.na(ena.out$depth) | ena.out$depth == ""),
  elevation = !(is.na(ena.out$elevation) | ena.out$elevation == "")
)
test_syntax <- data.frame(
  altitude = grepl("(0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?", ena.out$altitude) & test_pa$altitude,
  depth = grepl("(0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?", ena.out$depth) & test_pa$depth,
  elevation = grepl("[+-]?(0|((0\\.)|([1-9][0-9]*\\.?))[0-9]*)([Ee][+-]?[0-9]+)?", ena.out$elevation) & test_pa$elevation
)
ena.assess$z.coordinate <- ifelse(
  !apply(test_pa, 1, any),
  "no",
  ifelse(
    apply(test_syntax, 1, any),
    "syntax",
    "other"
  )
)
ena.assess$z.coordinate <- factor(ena.assess$z.coordinate, levels = c("syntax", "other", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + z.coordinate ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_z_coordinate_available.txt", sep = "\t", quote = F, row.names = F)


### check nominal length ####
sum(is.na(ena.out$nominal_length))
# 784692 NA
sum(ena.out$nominal_length == "", na.rm = T)
# no empty cells
quantile(ena.out$nominal_length, seq(0, 1, 0.05), na.rm = T)
# some very weird values (incl. some larger 10kbp and less than 50bp)
hist(ena.out$nominal_length[!is.na(ena.out$nominal_length) & ena.out$nominal_length < 2000], breaks = 1000)
abline(v = c(150, 250, 300, 500, 600), col = "red")
# suspicious peaks at typical read lengths
ena.assess$nominal.length <- ifelse(
  is.na(ena.out$nominal_length),
  "no",
  ifelse(
    ena.out$nominal_length %in% c(0, 150, 250, 300, 500, 600) | ena.out$nominal_length > 2000,
    "suspicious",
    "yes"
  )
)
ena.assess$nominal.length <- factor(ena.assess$nominal.length, levels = c("yes", "suspicious", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + nominal.length ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_nominal_length.txt", sep = "\t", quote = F, row.names = F)


### check target gene ####
sum(is.na(ena.out$target_gene))
sum(is.na(ena.out$library_construction_protocol))
# 0, 436 NA cells
sum(ena.out$target_gene == "")
sum(ena.out$library_construction_protocol == "", na.rm = T) # this will give 234 unique entries, many which have nothing to do with library construction
# 865051, 838821 empty cells
# no consistent terminology (no controlled vocabulary used) --> this field is not machine readable
ena.assess$target.gene <- ifelse(
  ena.out$target_gene != "", 
  "yes",
  ifelse(
    ena.out$sample_accession %in% accnos_xml_primers,
    "xml",
    "no"
  )
)
ena.assess$target.gene <- factor(ena.assess$target.gene, levels = c("yes", "xml", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + target.gene ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_target_gene.txt", sep = "\t", quote = F, row.names = F)


### check inconsistency in library preparation ####
# search was conducted for library source = METAGENOMIC and library strategy = AMPLICON
# there should not by any RT-PCR or random PCR library selection
table(ena.out$library_selection)
ena.assess$library.inconsistencies <- ifelse(ena.out$library_selection != "PCR", "library metadata inconsistent", "no inconsistencies")
ena.assess$library.inconsistencies <- factor(ena.assess$library.inconsistencies, levels = c("no inconsistencies", "library metadata inconsistent"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + library.inconsistencies ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_library_inconsistencies.txt", sep = "\t", quote = F, row.names = F)


### check correct use of ontology terms ####
sum(is.na(ena.out$environment_biome))
sum(is.na(ena.out$environment_material))
sum(is.na(ena.out$environment_feature))
sum(is.na(ena.out$broad_scale_environmental_context))
sum(is.na(ena.out$local_environmental_context))
sum(is.na(ena.out$environmental_medium))
# no NA cells
sum(ena.out$environment_biome == "")
sum(ena.out$environment_material == "")
sum(ena.out$environment_feature == "")
sum(ena.out$broad_scale_environmental_context == "")
sum(ena.out$local_environmental_context == "")
sum(ena.out$environmental_medium == "")
all.equal(ena.out$environment_biome, ena.out$broad_scale_environmental_context)
all.equal(ena.out$environment_feature, ena.out$local_environmental_context)
all.equal(ena.out$environment_material, ena.out$environmental_medium)
# content the same for old and new MIxS parameter names
# stick to new names

# values should have the following format:
#   name [ENVO:ID]|name [ENVO:ID]

# download obo of latest stable ENVO release from:
# https://github.com/EnvironmentOntology/envo
# wget https://raw.githubusercontent.com/EnvironmentOntology/envo/master/envo.obo
# extract number and label
#   grep -A2 "^\[Term\]" envo.obo | sed '/^--$/d' | grep -v "^\[Term\]" > envo_terms.txt
# this is a very crude way of parsing the obo file
# warning: there is an error in the obo file that GO:0016301 does not have a name (should be: kinase activity)
# manually adjust in envo_terms.txt file
envo <- data.frame(
  matrix(
    scan(what = "character", "envo/envo_terms.txt", sep = "\n"),
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
test <- ena.out %>% 
  select(broad_scale_environmental_context, local_environmental_context, environmental_medium) %>% 
  mutate(
    broad2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(broad_scale_environmental_context))))),
    local2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(local_environmental_context))))),
    medium2 = gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(environmental_medium)))))
  )
# run pattern match in bash with grep (faster)
write.table(envo, "envo/envo_formatted.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(test, "envo/ena_out_envo.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# bash code (using the formatted ENVO names written to file above)
# matches to ENVO ID
# cut -f1 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > broad_syntax.txt
# cut -f2 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > local_syntax.txt
# cut -f3 ena_out_envo.txt | grep -F -n -a -f <(cut -f3 envo_formatted.txt | sed -e 's/^/\[/' -e 's/$/\]/' | paste -d' ' <(cut -f5 envo_formatted.txt) -) > medium_syntax.txt

# cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > broad_id.txt
# cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > local_id.txt
# cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > medium_id.txt

# matches to ENVO name
# cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > broad_name.txt
# cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > local_name.txt
# cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > medium_name.txt

# matches to simplified ENVO name
# cut -f4 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > broad_name2.txt
# cut -f5 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > local_name2.txt
# cut -f6 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > medium_name2.txt

# read grep output
tmp.names <- matrix(
  c("envo/broad_syntax.txt", "envo/broad_id.txt", "envo/broad_name.txt", "envo/broad_name2.txt",
    "envo/local_syntax.txt", "envo/local_id.txt", "envo/local_name.txt", "envo/local_name2.txt",
    "envo/medium_syntax.txt", "envo/medium_id.txt", "envo/medium_name.txt", "envo/medium_name2.txt"),
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
ena.assess$env_broad <- ifelse(ena.out$broad_scale_environmental_context == "", "no", "filled")
ena.assess$env_broad[ena.in.envo$broad$name2$rn] <- "approx"
ena.assess$env_broad[unique(c(ena.in.envo$broad$name$rn, ena.in.envo$broad$id$rn))] <- "name_id"
ena.assess$env_broad[ena.in.envo$broad$syntax$rn] <- "syntax"
ena.assess$env_broad <- factor(ena.assess$env_broad, levels = c("syntax", "name_id", "approx", "filled", "no"))

ena.assess$env_local <- ifelse(ena.out$local_environmental_context == "", "no", "filled")
ena.assess$env_local[ena.in.envo$local$name2$rn] <- "approx"
ena.assess$env_local[unique(c(ena.in.envo$local$name$rn, ena.in.envo$local$id$rn))] <- "name_id"
ena.assess$env_local[ena.in.envo$local$syntax$rn] <- "syntax"
ena.assess$env_local <- factor(ena.assess$env_local, levels = c("syntax", "name_id", "approx", "filled", "no"))

ena.assess$env_medium <- ifelse(ena.out$environmental_medium == "", "no", "filled")
ena.assess$env_medium[ena.in.envo$medium$name2$rn] <- "approx"
ena.assess$env_medium[unique(c(ena.in.envo$medium$name$rn, ena.in.envo$medium$id$rn))] <- "name_id"
ena.assess$env_medium[ena.in.envo$medium$syntax$rn] <- "syntax"
ena.assess$env_medium <- factor(ena.assess$env_medium, levels = c("syntax", "name_id", "approx", "filled", "no"))

tmp <- cast(
  ena.assess,
  "MIxS + gfbio + env_broad ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_env_broad.txt", sep = "\t", quote = F, row.names = F)

tmp <- cast(
  ena.assess,
  "MIxS + gfbio + env_local ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_env_local.txt", sep = "\t", quote = F, row.names = F)

tmp <- cast(
  ena.assess,
  "MIxS + gfbio + env_medium ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_env_medium.txt", sep = "\t", quote = F, row.names = F)


### check status environmental sample ####
table(ena.out$environmental_sample, useNA = "ifany")
# before, there were no NA calls
# now, there are 880807 NA cells
# I assume that the entries for this parameter were revised by ENA
ena.out %>% filter(!environmental_sample) %>% View() # mostly negative controls
ena.out %>% filter(environmental_sample) %>% View() # should be mostly env, but also some artificial seawater samples included
# check definition of parameter
# at the moment, parameter is as useless for data mining as before (although less false information included)
# maybe this parameter can be fixed with predictive AI methods?

ena.assess$env_sample <- ifelse(
  is.na(ena.out$environmental_sample),
  "not available",
  ifelse(
    ena.out$environmental_sample,
    "yes",
    "no"
  )
)
ena.assess$env_sample <- factor(ena.assess$env_sample, levels = c("yes", "no", "not available"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + env_sample ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_env_sample.txt", sep = "\t", quote = F, row.names = F)


### check publication available ####
ena.assess$publication.available <- ifelse(ena.out$study_accession %in% c(xref.out$V1, xref.out$V2), "xref found", "not found")
ena.assess$publication.available <- factor(ena.assess$publication.available, levels = c("xref found", "not found"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + publication.available ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_publication_available.txt", sep = "\t", quote = F, row.names = F)


### visualization ####

pdf("summaries/summary_all_runs.pdf", width = 7, height = 6, onefile = T)
for(i in 7:18) {
  par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), oma = c(0, 0, 4, 0))
  tmp_all <- cast(
    ena.assess,
    paste0(colnames(ena.assess)[i], " ~ year.created"), 
    value = "study.accnos",
    fun.aggregate = "length",
    fill = 0,
    add.missing = T
  )
  barplot(
    as.matrix(
      data.frame(tmp_all[, -1], check.names = F)
    ),
    las = 2,
    legend.text = tmp_all[, 1],
    args.legend = list(x = "topleft", bty = "n")
  )
  barplot(
    prop.table(
      as.matrix(
        data.frame(tmp_all[, -1], check.names = F)
      ),
      2
    ) * 100,
    las = 2
  )
  tmp_mixs <- cast(
    ena.assess[ena.assess$MIxS, ],
    paste0(colnames(ena.assess)[i], " ~ year.created"), 
    value = "study.accnos",
    fun.aggregate = "length",
    fill = 0,
    add.missing = T
  )
  barplot(
    as.matrix(
      data.frame(tmp_mixs[, -1], check.names = F)
    ),
    las = 2
  )
  barplot(
    prop.table(
      as.matrix(
        data.frame(tmp_mixs[, -1], check.names = F)
      ),
      2
    ) * 100,
    las = 2
  )
  title(main = colnames(ena.assess)[i], outer = T, line = 2)
}
dev.off()


head(ena.assess)
ena.assess$tmp <- ena.assess$no.data == "no_data" | 
  ena.assess$collection.date.available == "no" |
  ena.assess$lat_lon == "no" | 
  ena.assess$target.gene == "no" | 
  ena.assess$library.inconsistencies == "library metadata inconsistent" |
  ena.assess$env_broad %in% c("filled", "no") |
  ena.assess$env_local %in% c("filled", "no") |
  ena.assess$env_medium %in% c("filled", "no")

tmp_all <- cast(
  ena.assess,
  "tmp ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
barplot(
  as.matrix(
    data.frame(tmp_all[, -1], check.names = F)
  ),
  las = 2,
  legend.text = tmp_all[, 1],
  args.legend = list(x = "topleft", bty = "n")
)

# for SAME poster ####
same_colors <- list(
  lat_lon = c("blue", "gold", "darkred"),
  collection.date.available = c("blue", "dodgerblue", "gold", "darkred"),
  target.gene = c("blue", "gold", "darkred"),
  env_broad = c("blue", "dodgerblue", "gold", "red", "darkred"),
  env_local = c("blue", "dodgerblue", "gold", "red", "darkred"),
  env_medium = c("blue", "dodgerblue", "gold", "red", "darkred"),
  publication.available = c("blue", "gold"),
  env_sample = c("blue", "dodgerblue", "darkred")
)
options(scipen=999)
for(i in names(same_colors)) {
  tmp_all <- cast(
    ena.assess[!ena.assess$MIxS, ],
    paste0(i, " ~ year.created"), 
    value = "study.accnos",
    fun.aggregate = "length",
    fill = 0,
    add.missing = T
  )
  tmp_mixs <- cast(
    ena.assess[ena.assess$MIxS, ],
    paste0(i, " ~ year.created"), 
    value = "study.accnos",
    fun.aggregate = "length",
    fill = 0,
    add.missing = T
  )
  png(paste0("summaries/same17_", i, "_v2.png"), width = 12, height = 5, units = "in", res = 600)
  par(mfrow = c(1, 2), mar = c(6, 5, 0.5, 0.5), ann = F)
  plot_ena_summaries(
    dat_all = tmp_all,
    dat_mixs = tmp_mixs,
    col = same_colors[[i]],
    ymax_all = 100000,
    ymax_mixs = 80000
  )
  dev.off()
}

tmp <- cast(
  ena.assess,
  "MIxS ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
