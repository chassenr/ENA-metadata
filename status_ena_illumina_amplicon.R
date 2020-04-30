# selecting samples for oligotyping Sulfurimonas
setwd("D:/Zmt_cloud/Data_mining")
require(data.table)
require(reshape)
require(rentrez)
require(XML)
# save.image("ENA_output/status_ena_illumina_amplicon.Rdata")
load("ENA_output/status_ena_illumina_amplicon.Rdata")

# OPEN QUESTION?
# take all amplicons studies, only ecological metagenomes, all metagenomes but exclude human?
# decide this first before looking into pubmed
# analysis at run or study level? I would prefer run level.
# Showcase GFBio? Only very few runs submitted this way compared to rest
# Is controled vocabulary for biome, material, feature mandatory?


# read new search results
ena.out <- fread(
  "ENA_output/20200428_results_read_run_tsv.txt",
  h = T,
  sep = "\t",
  quote = ""
)

# which columns are completely empty
empty.metadata <- colnames(ena.out)[apply(ena.out, 2, function(x) sum(is.na(x) | x == "") == length(x))]

# remove empty columns
ena.out <- ena.out[, (empty.metadata):=NULL]
dim(ena.out)
# 985311 entries, 111 fields with data

# retrieve taxid lineage
sum(is.na(ena.out$tax_id))
sum(ena.out$tax_id == "")
# taxid is always provided
taxid.uniq <- data.frame(
  uid = unique(ena.out$tax_id)
)
# retrieve lineage for each taxid
taxid.uniq$path <- c(NA)
taxid.uniq <- taxid.uniq[, c("path", "uid")]
for (i in 753:nrow(taxid.uniq)) {
  print(i)
  tmp <- scan(
    text = entrez_fetch("taxonomy", id = taxid.uniq$uid[i], rettype = "full", retmode = "xml", parsed=TRUE),
    what = "character",
    sep = "\n"
  )
  taxid.uniq[i, "path"] <- gsub("</Lineage>$", "", gsub("^ *<Lineage>", "", grep("<Lineage>", tmp, value = T)))
}
# classify into:
#   ecological metagenomes
#   organismal metagenomes
#   environmental sample?
#   other
taxid.uniq$category <- c()
taxid.uniq$category[grepl("ecological metagenomes", taxid.uniq$path)] <- "ecological metagenomes"
taxid.uniq$category[grepl("organismal metagenomes", taxid.uniq$path)] <- "organismal metagenomes"
taxid.uniq$category[grepl("environmental sample", taxid.uniq$path)] <- "environmental sample"
taxid.uniq$category[is.na(taxid.uniq$category)] <- "other"
# map to ena output
ena.out$taxid_cat <- taxid.uniq$category[match(ena.out$tax_id, taxid.uniq$uid)]
table(ena.out$taxid_cat)
# ecological metagenomes   environmental sample organismal metagenomes                  other 
# 325767                  15995                 509946                                  133603 


### target gene information available (run level) ####
sum(is.na(ena.out$target_gene))
table(ena.out$target_gene)
# no consistent terminology
sum(ena.out$target_gene == "")/nrow(ena.out) * 100
# 96.63% of runs no target gene specified
check.df <- data.frame(
  target.gene = ena.out$target_gene != "",
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "",
  taxid.cat = ena.out$taxid_cat
)
# comment: 
#   checklist information is not really useful because it only lists ERC[0-9]* and very few studies with MIxS,
#   although MIxS environmental package is provided for quite a few runs

# overview: target gene yes/no
check.summary <- cast(check.df, "target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -1]),
  col = c("red", "blue"),
  legend.text = check.summary$target.gene,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -1]), 2),
  col = c("red", "blue"),
  legend.text = check.summary$target.gene,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# apart from constantly increasing number of runs, no pattern over the years
# proportion of runs with target_gene info has been ~stable
# slight increase in studies with target_gene submitted between 2015 and 2019 
# (1.43% --> 3.40%)

# overview: target gene yes/no separated by taxid category
check.summary <- cast(check.df, "target.gene + taxid.cat ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in levels(check.summary$taxid.cat)) {
  barplot(
    as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
}
# no obviously different pattern between taxids: low proportion of runs with target_gene for organismal and ecological metagenomes
# slight increasing trend in TRUE proportion for ecological metagenomes from 2014 to 2019 (1.33 - 3.65%)
# overall similar proportion for organismal metagenomes, no trend
# only higher proportion of runs with target_gene which have 'environmental sample' in their taxonomic path
# however, compared to the number of runs for the other taxids, this is >1 order of magnitude less runs
# I would not split by taxid
# maybe focus only on ecological metagenomes or take all?

# overview: target gene yes/no separated by use of MIxS checklist
# all amplicons
check.summary <- cast(check.df, "target.gene + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}

# only ecological metagenomes
check.summary <- cast(check.df[check.df$taxid.cat == "ecological metagenomes", ], "target.gene + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# substantial increase of rund with target_gene info if checklists (MIxS environmental packages) are used
# this increase is more even if all studies are considered (not just ecological metagenomes)
# all runs 2015-2019: 1.8-20% with MIxS; ~1.8% no MIxS
# ecological metagenomes 2014-2019: 3.4-50% with MIxS; ~1% no MIxS
# however, only about 10% of runs use MIxS

# overview: target gene yes/no separated by use of GFBio
# all amplicons
check.summary <- cast(check.df, "target.gene + gfbio ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}

# only ecological metagenomes
table(check.df[check.df$gfbio, "taxid.cat"])
# as almost exclusively all GFBio runs are ecological metagenomes, above plot is sufficient
# too few GFBio submission compared to all runs
# about 50% of runs submitted via GFBio in 2017 and 2019 with target_gene infor
# I would not single out GFBio in the analysis
# rather I would mention in the conclusion that such brokerarge services can provide
# extremely valuable guidance for data submissions



### lat/lon available (run level) ####
sum(is.na(ena.out$lat))
sum(is.na(ena.out$lon))
sum(ena.out$lat == "", na.rm = T)
sum(ena.out$lon == "", na.rm = T)
sum(is.na(ena.out$lat) | is.na(ena.out$lon))/nrow(ena.out) * 100
# 36.14% of runs no lat/lon
check.df <- data.frame(
  lat.lon = !is.na(ena.out$lat) & !is.na(ena.out$lon),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "",
  taxid.cat = ena.out$taxid_cat
)

# overview: lat/lon yes/no
check.summary <- cast(check.df, "lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -1]),
  col = c("red", "blue"),
  legend.text = check.summary$lat.lon,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -1]), 2),
  col = c("red", "blue"),
  legend.text = check.summary$lat.lon,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# since 2014, majority of studies with lat/lon

# overview: lat/lon yes/no separated by taxid category
check.summary <- cast(check.df, "lat.lon + taxid.cat ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in levels(check.summary$taxid.cat)) {
  barplot(
    as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
}
# higher proportion of runs with lat/lon for ecological metagenomes
# this may also be related to privacy/confidentiality reasons
# should we exclude human samples from analysis?
# it may be easier to just focus on ecological metagenomes

# overview: lat/lon yes/no separated by use of MIxS checklist
# all amplicons
check.summary <- cast(check.df, "lat.lon + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}

# only ecological metagenomes
check.summary <- cast(check.df[check.df$taxid.cat == "ecological metagenomes", ], "lat.lon + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# for all amplicons, the use of checklists substantially increases the availability of lat/lon
# for ecological metagenomes, this difference is less pronounced
# in general, lat/lon mostly provided

# overview: lat/lon yes/no separated by use of GFBio
# all amplicons
check.summary <- cast(check.df, "lat.lon + gfbio ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# almost all GFBIO runs with lat/lon (exception in 2018)



### nominal length available (run level) ####
sum(is.na(ena.out$nominal_length))
sum(ena.out$nominal_length == "", na.rm = T)
sum(is.na(ena.out$nominal_length))/nrow(ena.out) * 100
# 78.32% of runs no nominal length, although mandatory for ENA, NCBI, DDBJ submissions
check.df <- data.frame(
  nom.len = !is.na(ena.out$nominal_length),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "",
  taxid.cat = ena.out$taxid_cat
)

# overview: nominal length yes/no
check.summary <- cast(check.df, "nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -1]),
  col = c("red", "blue"),
  legend.text = check.summary$nom.len,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -1]), 2),
  col = c("red", "blue"),
  legend.text = check.summary$nom.len,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# since 2013, constant decrease in proportion of studies with nominal length provided

# overview: nominal length yes/no separated by taxid category
check.summary <- cast(check.df, "nom.len + taxid.cat ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in levels(check.summary$taxid.cat)) {
  barplot(
    as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
}
# same pattern of increasing proportion without nominal length
# except for environmental sample in taxonomic path
# for those, nominal length almost always available
# but again, those only represent a small fraction of runs

# overview: nominal length yes/no separated by use of MIxS checklist
# all amplicons
check.summary <- cast(check.df, "nom.len + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}

# only ecological metagenomes
check.summary <- cast(check.df[check.df$taxid.cat == "ecological metagenomes", ], "nom.len + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# for all amplicons, the use of checklists substantially increases the availability of nominal length
# regardless of whether all or only ecological mategenomes are considered
# the difference is between 90% no nominal length (without MIxS) vs 10% (with MIxS)

# overview: nominal length yes/no separated by use of GFBio
# all amplicons
check.summary <- cast(check.df, "nom.len + gfbio ~ year.created", value = "study.accnos", fun.aggregate = "length")
# all GFBIO runs with nominal length

# for now, don't check potentially wrong entries
# just mention that there is some confusion regarding definition of nominal length
# check for subset of studies for Sulfurimonas oligotyping (easier to confirm)



### correct use of ontology terms ####
# get advice from Pier?
# if provided, do these 3 fields have to contain controlled vocabulary?
sum(is.na(ena.out$environment_biome))
sum(is.na(ena.out$environment_material))
sum(is.na(ena.out$environment_feature))
sum(ena.out$environment_biome != "" & ena.out$environment_material != "" & ena.out$environment_feature != "")/nrow(ena.out) * 100
# 28.39% of runs with entries in all 3 fields

# download json of ENVO terms from:
# https://www.ebi.ac.uk/ols/ontologies/envo/terms
# extract number and label
# grep -A2 "^\[Term\]" envo.obo | sed '/^--$/d' | grep -v "^\[Term\]" > envo_terms.txt
envo <- data.frame(
  matrix(
    scan(what = "character", "envo_terms.txt", sep = "\n"),
    ncol = 2,
    byrow = T
  )
)
colnames(envo) <- c("id", "name")
sum(grepl("id\\:", envo$id)) == nrow(envo)
sum(grepl("name\\:", envo$name)) == nrow(envo)
envo[!grepl("name\\:", envo$name), ]
# remove the last of these entries
envo <- envo[-5945, ]
envo[!grepl("name\\:", envo$name), ]
# should be ok
envo$id2 <- gsub("^id\\: ", "", envo$id)
envo$name2 <- gsub("^name\\: ", "", envo$name)
# write to file (faster to check in bash?)
write.table(envo, "envo_formatted.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# in bash
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f28 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > biome.in.envo1
# sed -i 's/:/\t/' biome.in.envo1
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' | sed 's/:/_/' > tmp
# cut -f28 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > biome.in.envo2
# sed -i 's/:/\t/' biome.in.envo2
# cut -f4 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f28 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -nw -f tmp > biome.in.envo3
# sed -i 's/:/\t/' biome.in.envo3
# 
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f29 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > feature.in.envo1
# sed -i 's/:/\t/' feature.in.envo1
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' | sed 's/:/_/' > tmp
# cut -f29 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > feature.in.envo2
# sed -i 's/:/\t/' feature.in.envo2
# cut -f4 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f29 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -nw -f tmp > feature.in.envo3
# sed -i 's/:/\t/' feature.in.envo3
# 
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f30 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > material.in.envo1
# sed -i 's/:/\t/' material.in.envo1
# cut -f3 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' | sed 's/:/_/' > tmp
# cut -f30 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -n -f tmp > material.in.envo2
# sed -i 's/:/\t/' material.in.envo2
# cut -f4 ../envo_formatted.txt | tr '[:upper:]' '[:lower:]' > tmp
# cut -f30 20200428_results_read_run_tsv.txt | sed '1d' | tr '[:upper:]' '[:lower:]' | grep -nw -f tmp > material.in.envo3
# sed -i 's/:/\t/' material.in.envo3

# read grep output
temp.names <- c("biome.in.envo1", "biome.in.envo2", "biome.in.envo3")
biome.in.envo <- vector("list", length = 3)
for(i in 1:length(temp.names)) {
  biome.in.envo[[i]] <- read.table(
    paste("ENA_output/", temp.names[i], sep = ""),
    h = F,
    sep = "\t",
    stringsAsFactors = F
  )
}
biome.in.envo.rn <- sort(unique(unlist(sapply(biome.in.envo, function(x) x[, 1]))))
biome.in.envo.log <- vector("logical", length = nrow(ena.out))
biome.in.envo.log[biome.in.envo.rn] <- TRUE
biome.in.envo.log[!biome.in.envo.log] <- FALSE

temp.names <- c("feature.in.envo1", "feature.in.envo2", "feature.in.envo3")
feature.in.envo <- vector("list", length = 3)
for(i in 1:length(temp.names)) {
  feature.in.envo[[i]] <- read.table(
    paste("ENA_output/", temp.names[i], sep = ""),
    h = F,
    sep = "\t",
    stringsAsFactors = F
  )
}
feature.in.envo.rn <- sort(unique(unlist(sapply(feature.in.envo, function(x) x[, 1]))))
feature.in.envo.log <- vector("logical", length = nrow(ena.out))
feature.in.envo.log[feature.in.envo.rn] <- TRUE
feature.in.envo.log[!feature.in.envo.log] <- FALSE

temp.names <- c("material.in.envo1", "material.in.envo2", "material.in.envo3")
material.in.envo <- vector("list", length = 3)
for(i in 1:length(temp.names)) {
  material.in.envo[[i]] <- read.table(
    paste("ENA_output/", temp.names[i], sep = ""),
    h = F,
    sep = "\t",
    stringsAsFactors = F
  )
}
material.in.envo.rn <- sort(unique(unlist(sapply(material.in.envo, function(x) x[, 1]))))
material.in.envo.log <- vector("logical", length = nrow(ena.out))
material.in.envo.log[material.in.envo.rn] <- TRUE
material.in.envo.log[!material.in.envo.log] <- FALSE

# overall stats:
sum(biome.in.envo.log & material.in.envo.log & feature.in.envo.log)/nrow(ena.out) * 100
# 11.73% with ENVO compatible terms in all runs
sum(biome.in.envo.log & material.in.envo.log & feature.in.envo.log)/sum(ena.out$environment_biome != "" & ena.out$environment_material != "" & ena.out$environment_feature != "") * 100
# 41.31% with ENVO compatible terms in runs with entries in these 3 fields

# How to proceed? For now, comparison against all runs
check.df <- data.frame(
  envo = biome.in.envo.log & material.in.envo.log & feature.in.envo.log,
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "",
  taxid.cat = ena.out$taxid_cat
)

# overview: envo yes/no
check.summary <- cast(check.df, "envo ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -1]),
  col = c("red", "blue"),
  legend.text = check.summary$envo,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -1]), 2),
  col = c("red", "blue"),
  legend.text = check.summary$envo,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# used since 2014, but on average only small proportion of studies, decreasing over time

# overview: envo yes/no separated by taxid category
check.summary <- cast(check.df, "envo + taxid.cat ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in levels(check.summary$taxid.cat)) {
  barplot(
    as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
}
# slightly higher proportion for ecological metagenomes
# but decreasing trend in envo usage since 2015

# overview: envo yes/no separated by use of MIxS checklist
# all amplicons
check.summary <- cast(check.df, "envo + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}

# only ecological metagenomes
check.summary <- cast(check.df[check.df$taxid.cat == "ecological metagenomes", ], "envo + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# higher proportion of envo in runs with MIxS (~40-50%)
# the previously mentioned decreasing trend of envo usage is not observed in runs with MIxS

# overview: envo yes/no separated by use of GFBio
# all amplicons
check.summary <- cast(check.df, "envo + gfbio ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# almost all GFBio with correct envo terms
View(ena.out[ena.out$broker_name == "GFBIO" & !check.df$envo, ])



### declared as environmental sample ####
sum(is.na(ena.out$environmental_sample))
prop.table(table(ena.out$environmental_sample)) * 100
# only 1.3 of runs declared as environmental samples
# based on manually curated subset for Sulfurimonas I know that this is not correct

check.df <- data.frame(
  env.sample = ena.out$environmental_sample,
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "",
  taxid.cat = ena.out$taxid_cat
)

# overview: environmental sample yes/no
check.summary <- cast(check.df, "env.sample ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
barplot(
  as.matrix(check.summary[, -1]),
  col = c("red", "blue"),
  legend.text = check.summary$env.sample,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
barplot(
  prop.table(as.matrix(check.summary[, -1]), 2),
  col = c("red", "blue"),
  legend.text = check.summary$env.sample,
  args.legend = list(x = "topleft", cex = 0.7),
  names.arg = colnames(check.summary)[-1],
  las = 2,
  cex.axis = 0.7,
  cex.names = 0.7
)
# slighly higher proportion in 2017 - 2019

# overview: environmental sample yes/no separated by taxid category
check.summary <- cast(check.df, "env.sample + taxid.cat ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in levels(check.summary$taxid.cat)) {
  barplot(
    as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$taxid.cat == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = i
  )
}
# only high proportion in runs with environmental sample in taxonomic path

# overview: environmental sample yes/no separated by use of MIxS checklist
# all amplicons
check.summary <- cast(check.df, "env.sample + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}

# only ecological metagenomes
check.summary <- cast(check.df[check.df$taxid.cat == "ecological metagenomes", ], "env.sample + MIxS ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$MIxS)) {
  barplot(
    as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$MIxS == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("MIxS", i)
  )
}
# there are no ecological metagenomes declared as environmental samples
# otherwise, considering all samples --> no pattern

# overview: target gene yes/no separated by use of GFBio
# all amplicons
check.summary <- cast(check.df, "env.sample + gfbio ~ year.created", value = "study.accnos", fun.aggregate = "length")
par(mfrow = c(1, 2))
for(i in unique(check.summary$gfbio)) {
  barplot(
    as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
  barplot(
    prop.table(as.matrix(check.summary[check.summary$gfbio == i, -c(1, 2)]), 2),
    col = c("red", "blue"),
    legend.text = c("FALSE", "TRUE"),
    args.legend = list(x = "topleft", cex = 0.7),
    names.arg = colnames(check.summary)[-c(1, 2)],
    las = 2,
    cex.axis = 0.7,
    cex.names = 0.7,
    main = paste("gfbio", i)
  )
}
# there are no GFBIO submitted sample declared as environmental samples
# I know for sure that this is not correct



### retrieve pubmed ID ####
# maybe skip for large data set and only look at this for subset for Sulfurimonas


# write unique primary study accnos to file
write.table(
  unique(ena.out.1$study_accession),
  "ena_primary_study_accnos.txt",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)

# read primary study accnos without xref hits
xref.no.primary.study.accnos <- gsub(
  "_xref.txt",
  "",
  gsub(
    "./",
    "",
    read.table(
      "no_hits_primary.txt",
      sep = "\t",
      stringsAsFactors = F
    )[, 1],
    fixed = T
  ),
  fixed = T
)
# write unique secondary study accnos to file
write.table(
  unique(ena.out.1$secondary_study_accession[ena.out.1$study_accession %in% xref.no.primary.study.accnos]),
  "ena_secondary_study_accnos.txt",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)
# there seem to be more unique secondary study accessions than primary ones...
temp <- unique(ena.out.1[, c("study_accession", "secondary_study_accession")])
temp[temp$study_accession %in% names(which(table(temp$study_accession) > 1)), ]
# PRJDB2953 and PRJDB3826

# read secondary study accnos without xref hits
xref.no.secondary.study.accnos <- gsub(
  "_xref.txt",
  "",
  gsub(
    "./",
    "",
    read.table(
      "no_hits_secondary.txt",
      sep = "\t",
      stringsAsFactors = F
    )[, 1],
    fixed = T
  ),
  fixed = T
)

# summary on study level
length(unique(ena.out.1$study_accession[ena.out.1$secondary_study_accession %in% xref.no.secondary.study.accnos]))/length(unique(ena.out.1$study_accession)) * 100
# 82.97573% of studies no retrievable pubmed ID










### subset of studies for Sulfurimonas oligotypting ####
oligotyping.studies <- read.table("ena_subset_sulfurimonas.txt", h = T, sep = "\t", stringsAsFactors = F)
ena.sub <- ena.out.1[ena.out.1$study_accession %in% oligotyping.studies$study_accession, ]
# which one are not included:
oligotyping.studies[!oligotyping.studies$study_accession %in% ena.out.1$study_accession, ]
# because of wrong library strategy: PRJEB14127, PRJEB11384
# because of wrong library layout: PRJEB15554, PRJEB18774, PRJNA299110, PRJNA282077, PRJEB20733
# because of wrong library selection: PRJEB10576

# were reads archived correctly?
sum(oligotyping.studies$Demultiplexed & oligotyping.studies$primer.clipped & oligotyping.studies$not.merged)
# 8 studies correct
table(oligotyping.studies$Demultiplexed)
# 1 study not demultiplexed
table(oligotyping.studies$primer.clipped)
# 28 studies not primer clipped
table(oligotyping.studies$not.merged)
# 8 studies already merged

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
# settelment experiments in natural locations also classified as environmental sample
sum(oligotyping.studies$env.sample.ena == oligotyping.studies$env.sample.real)
# 3 studies with correct value for environmental sample

# can publication be linked if available?
oligotyping.studies.noxref <- oligotyping.studies[oligotyping.studies$study_accession %in% ena.out.1$study_accession[ena.out.1$secondary_study_accession %in% xref.no.secondary.study.accnos], ]
nrow(oligotyping.studies.noxref)
# 19 studies without xref hits
sum(oligotyping.studies.noxref$publication.available)
# for at least 13 of those publication is actually available

###############################



