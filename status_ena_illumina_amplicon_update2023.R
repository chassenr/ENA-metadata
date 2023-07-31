# Assessing status of metadata on ENA for Illumina amplicon studies 
# from metagenomic material (taxid: ecological metagenomes)

### prepare environment and read data ####

setwd("C:/Users/chassenrueck/Documents/Additional_projects/Data_mining_commentary/Update_2023/")
require(data.table)
require(reshape)
require(XML)
require(dplyr)
require(tidyverse)
require(janitor)
require(furrr)

# save.image("status_ena_illumina_amplicon_update2023.Rdata")
load("status_ena_illumina_amplicon_update2023.Rdata")

# set multithreading options
# plan(multicore, workers = 48)

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
ena.assess$collection.date.available <- ifelse(
  !(is.na(ena.out$collection_date) | ena.out$collection_date == ""),
  "yes",
  ifelse(
    !(ena.out$collection_date_start == "" | ena.out$collection_date_end == ""),
    "period",
    "no"
  )
)
ena.assess$collection.date.available <- factor(ena.assess$collection.date.available, levels = c("yes", "period", "no"))
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
ena.assess$lat_lon <- ifelse(!is.na(ena.out$lat) & !is.na(ena.out$lon), "yes", "no")
ena.assess$lat_lon <- factor(ena.assess$lat_lon, levels = c("yes", "no"))
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
  ena.out$country == "",
  "no",
  ifelse(
    apply(test, 1, any),
    "insdc_country",
    "other"
  )
)
ena.assess$country.available <- factor(ena.assess$country.available, levels = c("insdc_country", "other", "no"))
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


### check z-ccordinates ####
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
    ena.out$nominal_length %in% c(0, 150, 250, 300, 500, 600),
    "question",
    "yes"
  )
)
ena.assess$nominal.length <- factor(ena.assess$nominal.length, levels = c("yes", "question", "no"))
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
ena.assess$target.gene <- ifelse(ena.out$target_gene != "", "yes", "no")
ena.assess$target.gene <- factor(ena.assess$target.gene, levels = c("yes", "no"))
tmp <- cast(
  ena.assess,
  "MIxS + gfbio + target.gene ~ year.created", 
  value = "study.accnos",
  fun.aggregate = "length",
  fill = 0,
  add.missing = T
)
write.table(tmp, "summaries/summary_target_gene.txt", sep = "\t", quote = F, row.names = F)
ena.assess %>% filter(gfbio & MIxS & !target.gene & year.created == 2023) %>% View()
# in some cases, target subfragment and pcr primers are provided, but not target gene,
# since the other 2 parameters are not indexed, the information is not accessible in the TSV output


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
