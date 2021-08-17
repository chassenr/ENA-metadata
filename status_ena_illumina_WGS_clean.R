# Assessing status of metadata on ENA for Illumina WGS studies 
# from metagenomic material (taxid: ecological metagenomes)

### prepare environment and read data ####

setwd("C:/Users/tobia/Desktop/ENA Paper/")
require(data.table)
require(reshape)
require(XML)
require(dplyr)
require(tidyverse)
require(janitor)
require(venn)

# save.image("ENA_out/status_ena_illumina_WGS.Rdata")
load("ENA_out/status_ena_illumina_WGS.Rdata")

# read ENA search results
ena.out <- fread(
  "results_read_run_tsv_WGS_metagenomes.txt",
  h = T,
  sep = "\t",
  quote = ""
)

# which columns are completely empty
empty.metadata <- colnames(ena.out)[apply(ena.out, 2, function(x) sum(is.na(x) | x == "") == length(x))]

# remove empty columns
ena.out <- ena.out[, (empty.metadata):=NULL]
dim(ena.out)
# 29253 entries, 95 fields with data

# check taxid
sum(is.na(ena.out$tax_id))
sum(ena.out$tax_id == "")
# taxid is always provided

# retrieving additional data from the sample xml records
# sample XML were directly downloaded from the results interface of the ENA search

sampleXML <- xmlParse("results_read_run_xml_WGS_metagenomes_20210108.xml")
sampleXML.top <- xmlRoot(sampleXML)
sampleXML.list <- map(
  1:xmlSize(sampleXML.top),
  function(X) {
    xmlToDataFrame(sampleXML.top[[X]][['SAMPLE_ATTRIBUTES']], stringsAsFactors = FALSE,) %>% 
      mutate_all(~type.convert(., as.is = T))
  }
)
sampleXML.accnos <- map_chr(1:xmlSize(sampleXML.top), function(X) xmlValue(sampleXML.top[[X]][['IDENTIFIERS']][['PRIMARY_ID']]))
sampleXML.list.df <- vector("list", length = length(sampleXML.list))
for(i in 1:length(sampleXML.list.df)) {
  suppressMessages(
    sampleXML.list.df[[i]] <- t(sampleXML.list[[i]][, 1:2]) %>% 
      row_to_names(row_number = 1) %>% 
      as_tibble(.name_repair = "universal") %>% 
      clean_names()
  )
}
names(sampleXML.list.df) <- sampleXML.accnos

# extract all potential column names that may indicate that runs were not WGS, but AMPLICON
sampleXML.colnames <- sort(unique(unlist(lapply(sampleXML.list.df, colnames))))
write.table(sampleXML.colnames, "sample_xml_metadata.txt", sep = "\t", row.names = F, col.names = F)
# potentially useful columns
check.columns <- c(
  "amplicon_library_type",
  "pcr_primers" 
)

# parse matrix for these columns
sampleXML.parsed <- matrix(NA, nrow = length(sampleXML.list.df), ncol = length(check.columns))
rownames(sampleXML.parsed) <- names(sampleXML.list.df)
colnames(sampleXML.parsed) <- check.columns
for(i in 1:nrow(sampleXML.parsed)) {
  if(any(colnames(sampleXML.list.df[[i]]) %in% check.columns)) {
    tmp.select <- colnames(sampleXML.list.df[[i]])[colnames(sampleXML.list.df[[i]]) %in% check.columns]
    sampleXML.parsed[i, tmp.select] <- as.matrix(sampleXML.list.df[[i]])[1, tmp.select] 
  }
}
# does not look like primer information
# however, pyrosequencing (+Solexa) was mentioned... 

# also inspect fields related to lat/lon
check.columns <- c(
  "geographic_location_latitude",
  "geographic_location_longitude",
  "lat_lon",
  "latitude_and_longitude",
  "latitude_end",
  "latitude_start",
  "latlon",
  "longitude_end",
  "longitude_start"
)
sampleXML.parsed2 <- matrix(NA, nrow = length(sampleXML.list.df), ncol = length(check.columns))
rownames(sampleXML.parsed2) <- names(sampleXML.list.df)
colnames(sampleXML.parsed2) <- check.columns
for(i in 1:nrow(sampleXML.parsed2)) {
  if(any(colnames(sampleXML.list.df[[i]]) %in% check.columns)) {
    tmp.select <- colnames(sampleXML.list.df[[i]])[colnames(sampleXML.list.df[[i]]) %in% check.columns]
    sampleXML.parsed2[i, tmp.select] <- as.matrix(sampleXML.list.df[[i]])[1, tmp.select] 
  }
}
# sampleXML.parsed2.summary <- apply(sampleXML.parsed2, 2, table, useNA = "always")
all.equal(ena.out$secondary_sample_accession, rownames(sampleXML.parsed2))
ena.out$add_info_lat_lon <- apply(sampleXML.parsed2, 1, function(x) any(!is.na(x)))


# as further information about the library preparation may be included in the experiment XML
# retrieve experiment XML for all cases

# bash code (using the TSV ENA search output to extract experiment accession numbers):
#   cut -f1,33,77 results_read_run_tsv_WGS_metagenomes.txt > results_read_run_tsv_sample_exp_run.txt
#   sed '1d' results_read_run_tsv_sample_exp_run.txt | cut -f2 | sed 's/^/https\:\/\/www\.ebi\.ac\.uk\/ena\/browser\/api\/xml\//' | parallel -j200 'wget -qO- {}' | sed '/<\/EXPERIMENT_SET>/,/<EXPERIMENT_SET>/d' >> experiment.xml
#   cp experiment.xml experiment_ori.xml
#   echo "</EXPERIMENT_SET>" >> experiment.xml

expXML <- xmlParse("ENA_out/experiment.xml")
expXML.top <- xmlRoot(expXML)
expXML.list <- c()
for(i in 1:xmlSize(expXML.top)) {
  if(is.null(expXML.top[[i]][['DESIGN']][['LIBRARY_DESCRIPTOR']][['LIBRARY_CONSTRUCTION_PROTOCOL']])) {
    expXML.list[i] <- NA
  } else {
    tmp <- xmlToList(expXML.top[[i]][['DESIGN']][['LIBRARY_DESCRIPTOR']][['LIBRARY_CONSTRUCTION_PROTOCOL']])
    expXML.list[i] <- ifelse(is.null(tmp), NA, tmp) 
  }
}
expXML.accnos <- map_chr(1:xmlSize(expXML.top), function(X) xmlValue(expXML.top[[X]][['IDENTIFIERS']][['PRIMARY_ID']]))
names(expXML.list) <- expXML.accnos
expXML.summary <- table(expXML.list, useNA = "always")
write.table(names(expXML.summary),  "experiment_xml_metadata.txt", sep = "\t", row.names = F, col.names = F)
# no indication that amplicon libraries were wrongly assigned to WGS

# define MIxS checklist names
mixs.checklists <- paste0("ERC0000", c(12:25, 31))

#####


### lat/lon available (run level) ####

# inspect data
sum(is.na(ena.out$lat))
sum(is.na(ena.out$lon))
# many NAs
sum(ena.out$lat == "", na.rm = T)
sum(ena.out$lon == "", na.rm = T)
# no empty cells
sum(!is.na(ena.out$lat) & !is.na(ena.out$lon))/nrow(ena.out) * 100
# 89.35152% of runs with lat/lon
sum((!is.na(ena.out$lat) & !is.na(ena.out$lon)) | ena.out$add_info_lat_lon)/nrow(ena.out) * 100
# 96.71829% with some kind of lat/lon info

# write summary for SI table
check.summary <- cast(
  check.df, 
  "MIxS + gfbio + lat.lon ~ year.created", 
  value = "study.accnos", 
  fun.aggregate = "length", 
  fill = 0, add.missing = T
)
write.table(check.summary, "summary_latlon.txt", sep = "\t", quote = F, row.names = F)

#####


### target gene information available (run level) ####

# this is to double-check that amplicon studies were not wrongly submitted as WGS
sum(is.na(ena.out$target_gene))
# no NA cells
table(ena.out$target_gene)
# many empty cells
# no consistent terminology (no controlled vocabulary used)
# this field is not machine readable
sum(ena.out$target_gene != "")/nrow(ena.out) * 100
# 0.02051072% of runs with target gene specified in indexed parameter
table(ena.out$target_gene)
# however those entries just say: Whole metagenome
# this is ok

#####


### nominal length available (run level mandatory parameter) ####

# inspect data
sum(is.na(ena.out$nominal_length))
# many NA
sum(ena.out$nominal_length == "", na.rm = T)
# no empty cells
sum(is.na(ena.out$nominal_length))/nrow(ena.out) * 100
# 67.58281% of runs no nominal length, although mandatory for ENA, NCBI, DDBJ submissions

# write summary for SI table
check.summary <- cast(
  check.df, 
  "MIxS + gfbio + nom.len ~ year.created", 
  value = "study.accnos", 
  fun.aggregate = "length", 
  fill = 0, 
  add.missing = T
)
write.table(check.summary, "summary_nomlen.txt", sep = "\t", quote = F, row.names = F)

#####


### correct use of ontology terms ####

# values for environment_biome, environment_material, and environment_feature
# should have the following format:
#   name [ENVO:ID]

# inspect data
sum(is.na(ena.out$environment_biome))
sum(is.na(ena.out$environment_material))
sum(is.na(ena.out$environment_feature))
# no NA cells
sum(ena.out$environment_biome != "" & ena.out$environment_material != "" & ena.out$environment_feature != "")/nrow(ena.out) * 100
# 34.48193% of runs with entries in all 3 fields

# download obo of latest stable ENVO release from:
# https://github.com/EnvironmentOntology/envo
# extract number and label:  
#   grep -A2 "^\[Term\]" envo.obo | sed '/^--$/d' | grep -v "^\[Term\]" > envo_terms.txt
envo <- data.frame(
  matrix(
    scan(what = "character", "ENVO/envo_terms.txt", sep = "\n"),
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
envo$name4 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(envo$name3)))))

# write to file to do grep in bash (faster)
# to 'rescue' values, 
#   all are converted to lower case, 
#   all non-alphanumeric characters are replaced with underscore
#   duplicate, trailing and leading underscores are removed
write.table(envo, "ENVO/envo_formatted.txt", sep = "\t", quote = F, row.names = F, col.names = F)
ena.out.envo <- ena.out[, c("environment_biome", "environment_material", "environment_feature")]
ena.out.envo$biome2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_biome)))))
ena.out.envo$material2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_material)))))
ena.out.envo$feature2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_feature)))))
write.table(ena.out.envo, "ena_out_envo.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# bash code (using the formatted ENVO names written to file above)
# matches to ENVO ID
#   cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > biome_id.txt
#   cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > material_id.txt
#   cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > feature_id.txt

# matches to ENVO name
#   cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > biome_name.txt
#   cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > material_name.txt
#   cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > feature_name.txt

# matches to simplified ENVO name
#   cut -f4 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > biome_name2.txt
#   cut -f5 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > material_name2.txt
#   cut -f6 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > feature_name2.txt

# read grep output
tmp.names <- matrix(
  c("ENVO/biome_id.txt", "ENVO/biome_name.txt", "ENVO/biome_name2.txt",
    "ENVO/material_id.txt", "ENVO/material_name.txt", "ENVO/material_name2.txt",
    "ENVO/feature_id.txt", "ENVO/feature_name.txt", "ENVO/feature_name2.txt"),
  nrow = 3,
  ncol = 3
)
ena.in.envo <- vector("list", length = 3)
names(ena.in.envo) <- c("biome", "material", "feature")
for(i in 1:length(ena.in.envo)) {
  ena.in.envo[[i]] <- vector("list", length = 3)
  names(ena.in.envo[[i]]) <- c("id", "name", "name2")
  for(j in 1:length(ena.in.envo[[i]])) {
    ena.in.envo[[i]][[j]] <- read.table(
      tmp.names[i, j],
      h = F,
      sep = "\t",
      stringsAsFactors = F
    )
    ena.in.envo[[i]][[j]]$rn <- as.numeric(gsub(":.*", "", ena.in.envo[[i]][[j]][, 1]))
  }
}

# parse output
check.df <- data.frame(
  biome = factor(ifelse(ena.out$environment_biome == "", "no", "filled"), levels = c("yes", "some", "filled", "no")),
  material = factor(ifelse(ena.out$environment_material == "", "no", "filled"), levels = c("yes", "some", "filled", "no")),
  feature = factor(ifelse(ena.out$environment_feature == "", "no", "filled"), levels = c("yes", "some", "filled", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "" & ena.out$checklist %in% mixs.checklists
)
for(i in names(ena.in.envo)) {
  check.df[unique(c(ena.in.envo[[i]]$name2$rn, ena.in.envo[[i]]$name$rn)), i] <- "some"
  check.df[ena.in.envo[[i]]$id$rn, i] <- "yes"
}
table(check.df$biome)
table(check.df$material)
table(check.df$feature)

# write summary for SI table
check.summary <- list(
  cast(check.df, "MIxS + gfbio + biome ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T),
  cast(check.df, "MIxS + gfbio + material ~ year.created", value = "study.accnos", fun.aggregate = "length",  fill = 0, add.missing = T),
  cast(check.df, "MIxS + gfbio + feature ~ year.created", value = "study.accnos", fun.aggregate = "length",  fill = 0, add.missing = T)
)
names(check.summary) <- names(ena.in.envo)
write.table(
  lapply(check.summary, function(x) { colnames(x)[3] <- "envo"; return(x) }) %>% do.call("rbind", .), 
  "summary_envo.txt",
  sep = "\t", 
  quote = F
)

#####


### declared as environmental sample ####

# inspect data
sum(is.na(ena.out$environmental_sample))
# no NA cells
prop.table(table(ena.out$environmental_sample)) * 100
# only 3.329573 of runs declared as environmental samples

check.df <- data.frame(
  env.sample = factor(ena.out$environmental_sample, levels = c("TRUE", "FALSE")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "" & ena.out$checklist %in% mixs.checklists
)

# overview: env sample yes/no separated by use of MIxS checklist
check.summary <- cast(check.df, "MIxS + env.sample ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$env.sample
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]),
    col = levels(droplevels(check.summary[check.summary$MIxS == i, "col"])),
    legend.text = levels(droplevels(check.summary$env.sample[check.summary$MIxS == i])),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(droplevels(check.summary[check.summary$MIxS == i, "col"])),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# over all years
check.all <- cast(check.df, "env.sample ~ MIxS", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
# none of the submissions using a MIxS environmental package, declared the sample as environmental!!!

#####

