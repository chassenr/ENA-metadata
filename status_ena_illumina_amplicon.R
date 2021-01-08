# Asessing status on metadata on ENA for Illumina amplicon studies from metagenomic material (taxid: ecological metagenomes)
setwd("C:/Users/chassenrueck/Documents/Additional_projects/Data_mining_commentary/Amplicon/")
require(data.table)
require(reshape)
# require(rentrez)
require(XML)
require(dplyr)
require(tidyverse)
# require(janitor)
require(venn)

# save.image("ENA_out/status_ena_illumina_amplicon.Rdata")
load("ENA_out/status_ena_illumina_amplicon.Rdata")

# read search results
ena.out <- fread(
  "ENA_out/results_read_run_tsv.txt",
  h = T,
  sep = "\t",
  quote = ""
)

# which columns are completely empty
empty.metadata <- colnames(ena.out)[apply(ena.out, 2, function(x) sum(is.na(x) | x == "") == length(x))]

# remove empty columns
ena.out <- ena.out[, (empty.metadata):=NULL]
dim(ena.out)
# 413849 entries, 104 fields with data

# check taxid
sum(is.na(ena.out$tax_id))
sum(ena.out$tax_id == "")
# taxid is always provided

# retrieving additional data from the xml records:
# sample ####
sampleXML <- xmlParse("ENA_out/ena_read_run_20201213-1900.xml")
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
for(i in 173004:length(test)) {
  suppressMessages(
    test[[i]] <- t(sampleXML.list[[i]][, 1:2]) %>% 
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
  "adapters",
  "amplicon",
  "amplicon_sequenced",
  "barcode_name",
  "forwar_primer_id",
  "forward_primer",
  "forward_primer_515f",
  "forward_primer_h279",
  "forward_primer_h280",
  "forward_primer_sequence_5_3",
  "forward_primer_ta_reuk454fwd1",
  "fwd_primer",
  "fwd_primer_name",
  "fwd_primer_seq",
  "gene",
  "gene_name",
  "gene_region",
  "gene_target",
  "linker_primer_sequence",
  "linkerprimer",
  "loci",
  "locus",
  "locus_target",
  "pcr_primers",
  "pcr_primers_2",
  "primer",
  "primer_f",
  "primer_f_reference",
  "primer_number",
  "primer_pair",
  "primer_r",
  "primer_r_reference",
  "primer_sequence_f",
  "primer_sequence_r",
  "primer_set",
  "primers",
  "primers_for_nif_h_amplification",
  "primers_target",
  "primerset",
  "rev_primer",
  "rev_primer_name",
  "rev_primer_seq",
  "reverse_primer",
  "reverse_primer_1612",
  "reverse_primer_1613",
  "reverse_primer_926r",
  "reverse_primer_id",
  "reverse_primer_sequence_5_3",
  "reverse_primer_ta_reuk_rev3",
  "rrna_region",
  "target",
  "target_gene",
  "target_gene_16",
  "target_gene_17",
  "target_region",
  "target_subfragment",
  "targeted",
  "targetgene"
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
# sampleXML.parsed.summary <- apply(sampleXML.parsed, 2, table, useNA = "always")
all.equal(ena.out$secondary_sample_accession, rownames(sampleXML.parsed))
ena.out$add_info_target_gene <- apply(sampleXML.parsed, 1, function(x) any(!is.na(x)))

# also check fields related to lat/lon
check.columns <- c(
  "estimated_latitude_of_seawater_source",
  "estimated_longitude_of_seawater_source",
  "geographic_location_latitude",
  "geographic_location_latitude_and_longitude",
  "geographic_location_longitude",
  "lat",
  "lat_lo_n",
  "lat_lon",
  "latitude",
  "latitude_and_longitude",
  "latitude_dd",
  "latitude_end",
  "latitude_start",
  "mixed_lat_lon",
  "sample_name_organism_host_isolation_source_collection_date_geo_loc_name_lat_lon",
  "soil_lat_lon_site_1",
  "soil_lat_lon_site_2",
  "longitude",
  "longitude_dd",
  "longitude_end",
  "longitude_start",
  "vv_lat",
  "vv_long"
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
check.exclude <- c(
  "1st PCR to amplify the target region and 2nd PCR to attach adaptors and index",
  "2-step PCR",
  "2 step PCR",
  "2 x 300 bp paired ends (PE) Illumina MiSeq v3",
  "2x300",
  "adapted from Kircher et al. (2011) and Gansauge & Meyer (2013)",
  "adapter ligation",
  "Amount of PCR products was estimated by gel image analysis, equal amounts (<100ng/sample) pooled and purified; equal amounts (measured fluorometrically) of subpools were pooled together and sequenced",
  "Barcoded amplicons",
  "Construction Protocol Details",
  "DNA fragmentation",
  "dual index",
  "Following the MiSeq library construction protocol",
  "Genomic DNA was randomly fragmented using nebulization, and a ~400-bp fraction (including adapters) was obtained by gel electrophoresis.",
  "illumina",
  "Illumina",
  "Illumina HiSeq 2500",
  "Illumina library",
  "Illumina MiSeq",
  "Illumina MiSeq 2x300",
  "Illumina MiSeq amplicon",
  "Illumina MiSeq Nextera XT and MiSeq Reagent Kit v3 (2*300)",
  "Illumina MiSeq reagent v3",
  "Illumina Nextera",
  "Illumina Nextera XT",
  "Illumina Nextera XT DNA Library Preparation Kit",
  "Illumina protocol",
  "Illumina TruSeq DNA library preparation protocol",
  "Illumina V2 kit                     2*250bp",
  "ILLUMINA_MISEQ",
  "Library construction protocol 1",
  "limited cycle amplicon",
  "Lysozyme-Achromopeptidase-ProteinaseK-phenol-chloroform method",
  "MiSeq",
  "MiSeq Reagent Kit v3",
  "MiSeq SOP",
  "MiSeq v2",
  "MiSeq V3 reagent kits",
  "missing",
  "na",
  "NA",
  "Nextera XT",
  "Nextera XT V3",
  "NextSeq_TruSeq RNA Library Prep Kit v2",
  "none provided",
  "null",
  "other",
  "Ovation Rapid DR Multiplex System 1-96 (NuGEN)",
  "Paired-end",
  "paired end",
  "PCR",
  "PCR based protocol",
  "PCR with illumina adapters",
  "PCR with Illumina adapters",
  "PCR with illumina dapters",
  "PCR with indexed primers",
  "PCR_with_Illumina_adapters",
  "samples mixed at equimolar concentration",
  "Standard Illumina protocol",
  "tailed PCR",
  "Tailed PCR",
  "The construction of the cDNA library was carried out with total RNA using the TruSeq RNA Sample Preparation Kit from Illumina",
  "The sequencing library is prepared by random fragmentation of the DNA or cDNA sample, followed",
  "This is a test upload",
  "TruDNA Seq library",
  "Truseq DNA PCR free Library Preparation Kit",
  "TruSeq PCR Free",
  "two-step PCR",
  "Two-Step PCR Approach",
  "Two-step tailed PCR",
  "two rounds of PCR, Dual-indexed libraries",
  "Two step PCR"
)
expXML.parsed <- expXML.list[!expXML.list %in% check.exclude & !is.na(expXML.list)]
sum(names(expXML.parsed) %in% ena.out$experiment_accession) == length(expXML.parsed)
ena.out$add_info_target_gene[ena.out$experiment_accession %in% names(expXML.parsed)] <- TRUE


### lat/lon available (run level) ####
sum(is.na(ena.out$lat))
sum(is.na(ena.out$lon))
# many NAs
sum(ena.out$lat == "", na.rm = T)
sum(ena.out$lon == "", na.rm = T)
# no empty cells
sum(!is.na(ena.out$lat) & !is.na(ena.out$lon))/nrow(ena.out) * 100
# 80.22177% of runs with lat/lon
sum((!is.na(ena.out$lat) & !is.na(ena.out$lon)) | ena.out$add_info_lat_lon)/nrow(ena.out) * 100
# 98.03141% with some kind of lat/lon info

check.df <- data.frame(
  lat.lon = factor(ifelse(!is.na(ena.out$lat) & !is.na(ena.out$lon), "yes", ifelse(ena.out$add_info_lat_lon, "some", "no")), levels = c("yes", "some", "no")),
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
# since 2014, majority of studies with (some kind of) lat/lon

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
# hardly any submissions without lat/lon information when MIxS is used
# BUT: some 15% not machine readable

# overview: lat/lon yes/no separated by use of GFBio
check.summary <- cast(check.df, "gfbio + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$lat.lon),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(check.summary$col),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# almost all GFBIO runs with good (machine readable) lat/lon
# over all years
check.all <- cast(check.df, "lat.lon ~ gfbio", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
View(ena.out[ena.out$study_accession %in% names(table(check.df[check.df$year.created == "2019" & check.df$gfbio & check.df$lat.lon == "no", "study.accnos"])), ])
# those are 2 negative controls

check.summary <- cast(check.df, "MIxS + gfbio + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_latlon.txt", sep = "\t", quote = F, row.names = F)

# summary:
#   after 2014 proportion without any lat/lon quickly decreased
#   since 2016 more of less constant proportion of runs without lat/lon over the years
#   hardly any submissions without lat/lon information when MIxS is used, but some 15% not machine readable
#   while only a few studies used a brokerage service (e.g.), this improved the availability of machine readable lat/lon


### target gene information available (run level) ####
sum(is.na(ena.out$target_gene))
# no NA cells
table(ena.out$target_gene)
# many empty cells
# no consistent terminology (no controlled vocabulary used)
# this field is not machine readable
sum(ena.out$target_gene != "")/nrow(ena.out) * 100
# 2.223516% of runs with target gene specified in indexed parameter
sum(ena.out$target_gene != "" | ena.out$add_info_target_gene)/nrow(ena.out) * 100
# 7.079635% of runs with some kind of info about target gene (hidden in sample or experiment XML)

check.df <- data.frame(
  target.gene = factor(ifelse(ena.out$target_gene != "", "yes", ifelse(ena.out$add_info_target_gene, "xml", "no")), levels = c("yes", "xml", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != ""
)
# comment: 
#   checklist information is not really useful because it only lists ERC[0-9]* and very few studies with MIxS,
#   although MIxS environmental package is provided for quite a few runs

# overview: target gene yes/no
check.summary <- cast(check.df, "target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$target.gene
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = check.summary$target.gene,
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
# most studies don't provide any kind of target gene information
# more than half that do, have this info hidden in xml
# no improvement over the years

# overview: target gene yes/no separated by use of MIxS checklist
check.summary <- cast(check.df, "MIxS + target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$target.gene
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$target.gene),
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
check.all <- cast(check.df, "target.gene ~ MIxS", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
# submissions with MIxS checklist (environmental package) have much higher availability of target gene info
# also more often in indexed parameter field

# overview: target gene yes/no separated by use of GFBio
check.summary <- cast(check.df, "gfbio + target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$target.gene
levels(check.summary$col) <- c("green", "yellow", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$target.gene),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(check.summary$col),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# no improvement compared to MIxS comparison
# except that almost exclusively indexed parameter field is used
check.all <- cast(check.df, "target.gene ~ gfbio", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)

check.summary <- cast(check.df, "MIxS + gfbio + target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_target_gene.txt", sep = "\t", quote = F, row.names = F)


### nominal length available (run level mandatory parameter) ####
sum(is.na(ena.out$nominal_length))
# many NA
sum(ena.out$nominal_length == "", na.rm = T)
# no empty cells
sum(is.na(ena.out$nominal_length))/nrow(ena.out) * 100
# 86.15075% of runs no nominal length, although mandatory for ENA, NCBI, DDBJ submissions
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
# number of studies providing this parameter has been constantly decreasing

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

# overview: nominal length yes/no separated by use of GFBio
check.summary <- cast(check.df, "gfbio + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
par(mfrow = c(2, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]),
    col = levels(check.summary$col),
    legend.text = levels(check.summary$nom.len),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2, ncol(check.summary))]), 2),
    col = levels(check.summary$col),
    names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# even better than MIxS
check.all <- cast(check.df, "nom.len ~ gfbio", value = "study.accnos", fun.aggregate = "length")
prop.table(as.matrix(check.all[, -1]), 2)
table(ena.out$nominal_length[ena.out$broker_name == "GFBIO"])

check.summary <- cast(check.df, "MIxS + gfbio + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
write.table(check.summary, "summary_nomlen.txt", sep = "\t", quote = F, row.names = F)


### correct use of ontology terms ####
# get advice from Pier?
# if provided, do these 3 fields have to contain controlled vocabulary?
# according to GFBio template, they should
sum(is.na(ena.out$environment_biome))
sum(is.na(ena.out$environment_material))
sum(is.na(ena.out$environment_feature))
# no NA cells
sum(ena.out$environment_biome != "" & ena.out$environment_material != "" & ena.out$environment_feature != "")/nrow(ena.out) * 100
# 25.11157% of runs with entries in all 3 fields

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
write.table(ena.out.envo, "ENVO/ena_out_envo.txt", sep = "\t", quote = F, row.names = F, col.names = F)

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
# decreasing trend of ENVO usage since 2015
# most of the time controlled vocabulary is used for material and feature,
# but hardly ever for biome

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
# almost all MIxS submission have entried for all 3 parameters, 
# but only about 80% (material, feature) and 12% (biome) using correct terminology

# overview: envo yes/no separated by use of GFBio
check.summary <- list(
  cast(check.df, "gfbio + biome ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "gfbio + material ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "gfbio + feature ~ year.created", value = "study.accnos", fun.aggregate = "length")
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 2]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "blue", "red")
}
par(mfrow = c(2, 2))
for(j in 1:length(check.summary)) {
  for(i in unique(check.summary[[j]]$gfbio)) {
    barplot(
      as.matrix(check.summary[[j]][check.summary[[j]]$gfbio == i, -c(1, 2, ncol(check.summary[[j]]))]),
      col = levels(droplevels(check.summary[[j]][check.summary[[j]]$gfbio == i, "col"])),
      legend.text = levels(droplevels(check.summary[[j]][check.summary[[j]]$gfbio == i, 2])),
      args.legend = list(x = "topleft", cex = 0.7),
      names.arg = colnames(check.summary[[j]])[-c(1, 2, ncol(check.summary[[j]]))],
      las = 2,
      cex.axis = 0.7,
      cex.names = 0.7,
      main = paste("gfbio", i)
    )
    barplot(
      prop.table(as.matrix(check.summary[[j]][check.summary[[j]]$gfbio == i, -c(1, 2, ncol(check.summary[[j]]))]), 2),
      col = levels(droplevels(check.summary[[j]][check.summary[[j]]$gfbio == i, "col"])),
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
  cast(check.df, "biome ~ gfbio", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "material ~ gfbio", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "feature ~ gfbio", value = "study.accnos", fun.aggregate = "length")
)
names(check.all) <- names(ena.in.envo)
lapply(check.all, function(x) prop.table(as.matrix(x[, -1]), 2))
# if provided, all submission through GFBio use correct ontology terms
# what happened in 2018
View(ena.out[ena.out$broker_name == "GFBIO" & ena.out$environment_biome == "", ])

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
# only 1.568446 of runs declared as environmental samples
# based on manually curated subset for Sulfurimonas I know that this is not correct

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


### subset of studies for Sulfurimonas oligotypting ####
oligotyping.studies <- read.table("ena_subset_sulfurimonas.txt", h = T, sep = "\t", stringsAsFactors = F)
ena.sub <- ena.out[ena.out$study_accession %in% oligotyping.studies$study_accession, ]
# which one are not included:
oligotyping.studies[!oligotyping.studies$study_accession %in% ena.out$study_accession, ]
# because of wrong library strategy: PRJEB14127, PRJEB11384
# because of wrong library layout: PRJEB15554, PRJEB18774, PRJNA299110, PRJNA282077, PRJEB20733
# because of wrong library selection: PRJEB10576
# Relevant studies were not included because of incorrectly provided metadata

# were reads archived correctly?
sum(oligotyping.studies$Demultiplexed & oligotyping.studies$primer.clipped & oligotyping.studies$not.merged)
# 8 studies correct
table(oligotyping.studies$Demultiplexed)
# 1 study not demultiplexed
table(oligotyping.studies$primer.clipped)
# 28 studies not primer clipped
table(oligotyping.studies$not.merged)
# 8 studies already merged
venn(oligotyping.studies[, 2:4])

# environmental sample declared accurately?
oligotyping.studies.runs <- read.table(
  "ena_subset_sulfurimonas_runs.txt",
  h = T,
  sep = "\t",
  stringsAsFactors = F,
  quote = "",
  comment.char = ""
)
table(oligotyping.studies.runs$environmental_sample)
# all declared as environmental sample FALSE
# this is not true
# manually checked
# study labeled as containing environmental samples, if at least one such samples was included
# other samples in the study may come from incubations etc.
# settlement experiments in natural locations also classified as environmental sample
sum(oligotyping.studies$env.sample.ena == oligotyping.studies$env.sample.real)
# 3 studies with correct value for environmental sample

# can publication be linked if available?
table(oligotyping.studies$publication.available)
# retrieve pubmed via XREF: see bash script
oligotyping.studies.noxref <- scan(what = "character", "ENA_xref/no_xref_study_accnos.txt")
# 19 studies without xref hits
sum(oligotyping.studies[oligotyping.studies$study_accession %in% oligotyping.studies.noxref, "publication.available"])
# for at least 14 of those, publication is actually available

###############################



