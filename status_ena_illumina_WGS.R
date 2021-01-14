# Asessing status on metadata on ENA for Illumina amplicon studies from metagenomic material (taxid: ecological metagenomes)
setwd("C:/Users/tobia/Desktop/ENA Paper/")

require(data.table)
require(reshape)
# require(rentrez)
require(XML)
require(dplyr)
require(tidyverse)
require(janitor)
require(venn)

#save.image("ENA_out/status_ena_illumina_WGS.Rdata")
load("ENA_out/status_ena_illumina_WGS.Rdata")

# read search results
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

# retrieving additional data from the xml records:
# sample ####
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

# extract all potential column names to search for useful primer information
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
#no specific inference of primer information
#but search query was explicit set for "library_selection"=RANDOM and not for parameter "random PCR" (which was found in xml output)
#and the search was explicitly restricted to "instrument_platform"=Illumina and pyrosequencing (+Solexa) was found


# sampleXML.parsed.summary <- apply(sampleXML.parsed, 2, table, useNA = "always")
#all.equal(ena.out$secondary_sample_accession, rownames(sampleXML.parsed))
#ena.out$add_info_target_gene <- apply(sampleXML.parsed, 1, function(x) any(!is.na(x)))

# also check fields related to lat/lon
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


# experiment ####
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


##same as target gene-xml view (=not specifically relevant)##


### lat/lon available (run level) ####
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

check.df <- data.frame(
  lat.lon = factor(ifelse(!is.na(ena.out$lat) & !is.na(ena.out$lon), "yes", ifelse(ena.out$add_info_lat_lon, "xml", "no")), levels = c("yes", "xml", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != ""
)

# overview: lat/lon yes/no
check.summary <- cast(check.df, "lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = check.summary$lat.lon,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -c(1, ncol(check.summary))]), 2),
  col = as.character(check.summary$col),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# since 2014, majority of studies with lat/lon info

# overview: lat/lon yes/no separated by use of MIxS checklist
check.summary <- cast(check.df, "MIxS + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$lat.lon),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(check.summary$col),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# over all years
check.all <- cast(check.df, "lat.lon ~ MIxS", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
# hardly any submissions with lat/lon information when MIxS is used


check.summary <- cast(check.df, "MIxS + gfbio + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_latlon.txt", sep = "\t", quote = F, row.names = F)

# summary:
#trend is pretty similiar to our example with amplicon studies, thus "non-prominent" result is shown in here 
#just some specific examples of wrong placement of information or rather wrongly provided information



### target gene information available (run level) ####
sum(is.na(ena.out$target_gene))
# no NA cells
table(ena.out$target_gene)
# many empty cells
# no consistent terminology (no controlled vocabulary used)
# this field is not machine readable
sum(ena.out$target_gene != "")/nrow(ena.out) * 100
# 0.02051072% of runs with target gene specified in indexed parameter

### nominal length available (run level mandatory parameter) ####
sum(is.na(ena.out$nominal_length))
# many NA
sum(ena.out$nominal_length == "", na.rm = T)
# no empty cells
sum(is.na(ena.out$nominal_length))/nrow(ena.out) * 100
# 67.58281% of runs no nominal length, although mandatory for ENA, NCBI, DDBJ submissions
check.df <- data.frame(
  nom.len = factor(ifelse(is.na(ena.out$nominal_length), "no", "yes"), levels = c("yes", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != ""
)

# overview: nominal length yes/no
check.summary <- cast(check.df, "nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = check.summary$nom.len,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -c(1, ncol(check.summary))]), 2),
  col = as.character(check.summary$col),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# number of studies providing this parameter has been nearly constant 

# overview: nominal length yes/no separated by use of MIxS checklist
check.summary <- cast(check.df, "MIxS + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$nom.len),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(check.summary$col),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# over all years
check.all <- cast(check.df, "nom.len ~ MIxS", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
# almost all MIxS studies with nominal length provided
# use of checklist definitely improved FAIRness

check.summary <- cast(check.df, "MIxS + gfbio + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_nomlen.txt", sep = "\t", quote = F, row.names = F)

#2019 check why filled studies are that representative
View(ena.out[check.df$material == "filled" & check.df$MIxS, ])
#https://www.ebi.ac.uk/ena/browser/view/PRJEB34634 this study makes it overrepresentative for 2019
#https://www.ebi.ac.uk/ena/browser/view/PRJEB31563 this one ist just a minor of it in 2019


#how many runs per WGS per study averaged
summary(c(table(ena.out$study_accession)))
#nominal length--> seems that read_length was used and not insert size
##----------------------------------------------------------------------------------------------
### correct use of ontology terms ####
# get advice from Pier?
# if provided, do these 3 fields have to contain controlled vocabulary?
# according to GFBio template, they should
sum(is.na(ena.out$environment_biome))
sum(is.na(ena.out$environment_material))
sum(is.na(ena.out$environment_feature))
# no NA cells
sum(ena.out$environment_biome != "" & ena.out$environment_material != "" & ena.out$environment_feature != "")/nrow(ena.out) * 100
# 34.48193% of runs with entries in all 3 fields

# download obo of latest stable ENVO release from:
# https://github.com/EnvironmentOntology/envo
# extract number and label
# grep -A2 "^\[Term\]" envo.obo | sed '/^--$/d' | grep -v "^\[Term\]" > envo_terms.txt
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
write.table(envo, "ENVO/envo_formatted.txt", sep = "\t", quote = F, row.names = F, col.names = F)
ena.out.envo <- ena.out[, c("environment_biome", "environment_material", "environment_feature")]
ena.out.envo$biome2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_biome)))))
ena.out.envo$material2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_material)))))
ena.out.envo$feature2 <- gsub("_$", "", gsub("^_", "", gsub("_{1,}", "_", gsub("[^[:alnum:]]", "_", tolower(ena.out.envo$environment_feature)))))
write.table(ena.out.envo, "ena_out_envo.txt", sep = "\t", quote = F, row.names = F, col.names = F)

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
  MIxS = ena.out$environmental_package != ""
)
for(i in names(ena.in.envo)) {
  check.df[unique(c(ena.in.envo[[i]]$name2$rn, ena.in.envo[[i]]$name$rn)), i] <- "some"
  check.df[ena.in.envo[[i]]$id$rn, i] <- "yes"
}
table(check.df$biome)
table(check.df$material)
table(check.df$feature)

# overview: envo yes/no
check.summary <- list(
  cast(check.df, "biome ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "material ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "feature ~ year.created", value = "study.accnos", fun.aggregate = "length")
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 1]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "blue", "red")
}
par(mfrow = c(3, 2))
for(i in 1:length(check.summary)) {
  barplot(
    as.matrix(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))]),
    col = as.character(check.summary[[i]]$col),
    legend.text = check.summary[[i]][, 1],
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary[[i]])[-c(1, ncol(check.summary[[i]]))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = names(check.summary)[i]
  )
  barplot(
    prop.table(as.matrix(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))]), 2),
    col = as.character(check.summary[[i]]$col),
    names.arg = colnames(check.summary[[i]])[-c(1, ncol(check.summary[[i]]))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7
  )
}
# ENVO usage have not changed a lot



# overview: envo yes/no separated by use of MIxS checklist
check.summary <- list(
  cast(check.df, "MIxS + biome ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + material ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + feature ~ year.created", value = "study.accnos", fun.aggregate = "length")
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 2]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "orange", "blue", "red")
}
par(mfrow = c(2, 2))
for(j in 1:length(check.summary)) {
  for(i in unique(check.summary[[j]]$MIxS)) {
    barplot(
      as.matrix(check.summary[[j]][check.summary[[j]]$MIxS == i, -c(1, 2, ncol(check.summary[[j]]))]),
      col = levels(droplevels(check.summary[[j]][check.summary[[j]]$MIxS == i, "col"])),
      legend.text = levels(droplevels(check.summary[[j]][check.summary[[j]]$MIxS == i, 2])),
      args.legend = list(x = "topleft", cex = 0.7),
      names.arg = colnames(check.summary[[j]])[-c(1, 2, ncol(check.summary[[j]]))],
      las = 2,
      cex.axis = 0.7,
      cex.names = 0.7,
      main = paste("MIxS", i)
    )
    barplot(
      prop.table(as.matrix(check.summary[[j]][check.summary[[j]]$MIxS == i, -c(1, 2, ncol(check.summary[[j]]))]), 2),
      col = levels(droplevels(check.summary[[j]][check.summary[[j]]$MIxS == i, "col"])),
      names.arg = colnames(check.summary[[j]])[-c(1, 2, ncol(check.summary[[j]]))],
      las = 2,
      cex.axis = 0.7,
      cex.names = 0.7,
      main = names(check.summary)[j]
    )
  }
}
# over all years
check.all <- list(
  cast(check.df, "biome ~ MIxS", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "material ~ MIxS", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "feature ~ MIxS", value = "study.accnos", fun.aggregate = "length")
)
names(check.all) <- names(ena.in.envo)
lapply(check.all, function(x) prop.table(as.matrix(x[, -1]), 2))
# decreasing trend is not seen in MIxS compatible submissions (i.e. using MIxS environmental package)
# almost all MIxS submission have entried for all 3 parameters



check.summary <- list(
  cast(check.df, "MIxS + gfbio + biome ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + gfbio + material ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + gfbio + feature ~ year.created", value = "study.accnos", fun.aggregate = "length")
)
names(check.summary) <- names(ena.in.envo)
write.table(
  lapply(check.summary, function(x) { colnames(x)[3] <- "envo"; return(x) }) %>% do.call("rbind", .), 
  "summary_envo.txt",
  sep = "\t", 
  quote = F
)


### declared as environmental sample ####
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
  MIxS = ena.out$environmental_package != ""
)

# overview: lat/lon yes/no
check.summary <- cast(check.df, "env.sample ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$env.sample
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = check.summary$env.sample,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -c(1, ncol(check.summary))]), 2),
  col = as.character(check.summary$col),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# hardly any samples declared as environmental

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

# overview: env sample yes/no separated by use of GFBio
check.summary <- cast(check.df, "gfbio + env.sample ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$env.sample
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]),
    col = levels(droplevels(check.summary[check.summary$gfbio == i, "col"])),
    legend.text = levels(check.summary$env.sample),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(droplevels(check.summary[check.summary$gfbio == i, "col"])),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# over all years
check.all <- cast(check.df, "env.sample ~ gfbio", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
# same here (don't show additional to MIxS plot)

check.summary <- cast(check.df, "MIxS + gfbio + env.sample ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_env_sample.txt", sep = "\t", quote = F, row.names = F)


