# plotting status ENA metadata

### figure 1: latlon, target gene, nominal length ####

# object size 750x550

pdf("/Users/Tobi/Desktop/Tobi/ENA Paper/ENA_out/plot1.pdf", width = 7.5, height = 5.5) # this will be in inches

par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.5), oma = c(1.5, 2, 2, 0))

# plot 1:
check.df <- data.frame(
  lat.lon = factor(ifelse(!is.na(ena.out$lat) & !is.na(ena.out$lon), "yes", ifelse(ena.out$add_info_lat_lon, "xml", "no")), levels = c("yes", "xml", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != ""
)
check.summary <- cast(check.df, "lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = paste0(check.summary$lat.lon, " (", round(rowSums(check.summary[, -c(1, ncol(check.summary))])/sum(rowSums(check.summary[, -c(1, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 1,
  cex.names = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3
)
mtext("All runs",3, 1, outer = F)
mtext("Latitude/longitude",3, 1, outer = F, font = 2, side = 3, line = 0.5, at = 31.2)
#abline(h = 10000)
# plot 2:
check.summary <- cast(check.df, "MIxS + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
barplot(
  as.matrix(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))]),
  col = levels(check.summary$col),
  legend.text = paste0(levels(check.summary$lat.lon), " (", round(rowSums(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))])/sum(rowSums(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
  las = 2,
  cex.axis = 1,
  cex.names = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3
)
mtext("MIxS environmental package", 3, 1, outer = F)

# plot 3:
check.df <- data.frame(
  nom.len = factor(ifelse(is.na(ena.out$nominal_length), "no", "yes"), levels = c("yes", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != ""
)
check.summary <- cast(check.df, "nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = paste0(check.summary$nom.len, " (", round(rowSums(check.summary[, -c(1, ncol(check.summary))])/sum(rowSums(check.summary[, -c(1, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  las = 2,
  cex.axis = 1,
  cex.names = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3
)
mtext("Nominal Length",3, 1, outer = F, font = 2, side = 3, line = 0.8, at = 28.75)
#abline(h = 8000)
# plot 4:
check.summary <- cast(check.df, "MIxS + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length")
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
barplot(
  as.matrix(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))]),
  col = levels(check.summary$col),
  legend.text = paste0(levels(check.summary$nom.len), " (", round(rowSums(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))])/sum(rowSums(check.summary[check.summary$MIxS, -c(1, 2, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
  las = 2,
  cex.axis = 1,
  cex.names = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3
)
# labels
mtext("Submission year", 1, 0, outer = T)
mtext("Number of runs", 2, 0, outer = T)

dev.off()

### figure 2: envo terms ####

# object size 700x750
pdf("/Users/Tobi/Desktop/Tobi/ENA Paper/ENA_out/plot2.pdf", width = 7, height = 7.5) # this will be in inches

par(mfcol = c(3, 2), mar = c(4, 4, 1, 0.5), oma = c(1.5, 2, 2, 0))

# plot 1-3:
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
for(i in 1:length(check.summary)) {
  barplot(
    as.matrix(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))]),
    col = as.character(check.summary[[i]]$col),
    legend.text = paste0(check.summary[[i]][, 1], " (", round(rowSums(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))])/sum(rowSums(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))])) * 100), "%)"),
    args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
    names.arg = colnames(check.summary[[i]])[-c(1, ncol(check.summary[[i]]))],
    las = 2,
    cex.axis = 1,
    cex.names = 1,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3
  )
  #abline(h = 8000)
  if(i == 1) {
    mtext("All runs", 3, 1, outer = F)
    mtext("Biome",3, 1, outer = F, font = 2, side = 3, line = 1.7, at = 20)
  }
  if(i == 2) {
       mtext("Material",3, 1, outer = F, font = 2, side = 3, line = 1.7, at = 21)
  }
  if(i == 3) {
    mtext("Feature",3, 1, outer = F, font = 2, side = 3, line = 1.7, at = 20.8)
  }
}

# plot 4-6:
check.summary <- list(
  cast(check.df, "MIxS + biome ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + material ~ year.created", value = "study.accnos", fun.aggregate = "length"),
  cast(check.df, "MIxS + feature ~ year.created", value = "study.accnos", fun.aggregate = "length")
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 2]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "blue", "red")
}
for(i in 1:length(check.summary)) {
  barplot(
    as.matrix(check.summary[[i]][check.summary[[i]]$MIxS, -c(1, 2, ncol(check.summary[[i]]))]),
    col = levels(droplevels(check.summary[[i]][check.summary[[i]]$MIxS, "col"])),
    legend.text = paste0(levels(droplevels(check.summary[[i]][check.summary[[i]]$MIxS, 2])), " (", round(rowSums(check.summary[[i]][check.summary[[i]]$MIxS, -c(1, 2, ncol(check.summary[[i]]))])/sum(rowSums(check.summary[[i]][check.summary[[i]]$MIxS, -c(1, 2, ncol(check.summary[[i]]))])) * 100), "%)"),
    args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
    names.arg = colnames(check.summary[[i]])[-c(1, 2, ncol(check.summary[[i]]))],
    las = 2,
    cex.axis = 1,
    cex.names = 1,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3
  )
  if(i == 1) {
    mtext("MIxS environmental package", 3, 1, outer = F)
    
  }
}
# labels
mtext("Submission year", 1, 0, outer = T)
mtext("Number of runs", 2, 0, outer = T)

dev.off()

