# plotting status ENA metadata

### figure 1: latlon, target gene, nominal length ####

pdf("ENA_fig_gps_gene_nomlen.pdf", width = 7, height = 7.5)
# object size 700x750
par(mfrow = c(3, 2), mar = c(4, 4, 1, 0.5), oma = c(1.5, 2, 2, 0), ann = F, yaxs = "i", lwd = 0.5)

# plot 1:
check.df <- data.frame(
  lat.lon = factor(ifelse(!is.na(ena.out$lat) & !is.na(ena.out$lon), "tsv", ifelse(ena.out$add_info_lat_lon, "xml", "no")), levels = c("tsv", "xml", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "" & ena.out$checklist %in% mixs.checklists
)
check.summary <- cast(check.df, "lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, max(table(check.df$year.created))),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
abline(h = 8000)
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = paste0(check.summary$lat.lon, " (", round(rowSums(check.summary[, -c(1, ncol(check.summary))])/sum(rowSums(check.summary[, -c(1, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
mtext("All cases", 3, 1, outer = F)
par(xpd = NA)
text(
  14.5, par("usr")[4],
  pos = 3,
  labels = "Geographic coordinates",
  font = 2,
  cex = 1.2
)
par(xpd = F)
# plot 2:
check.summary <- cast(check.df, "MIxS + lat.lon ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$lat.lon
levels(check.summary$col) <- c("green", "yellow", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, 8000),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
barplot(
  as.matrix(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))]),
  col = levels(check.summary$col),
  legend.text = paste0(
    levels(check.summary$lat.lon), 
    " (", 
    round(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])/sum(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])) * 100),
    "%)"
  ),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
mtext("MIxS environmental package", 3, 1, outer = F)
# plot 3:
check.df <- data.frame(
  target.gene = factor(ifelse(ena.out$target_gene != "", "tsv", ifelse(ena.out$add_info_target_gene, "xml", "no")), levels = c("tsv", "xml", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "" & ena.out$checklist %in% mixs.checklists
)
check.summary <- cast(check.df, "target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$target.gene
levels(check.summary$col) <- c("green", "yellow", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, max(table(check.df$year.created))),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
abline(h = 8000)
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = paste0(check.summary$target.gene, " (", round(rowSums(check.summary[, -c(1, ncol(check.summary))])/sum(rowSums(check.summary[, -c(1, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
par(xpd = NA)
text(
  14.5, par("usr")[4],
  pos = 3,
  labels = "Target gene",
  font = 2,
  cex = 1.2
)
par(xpd = F)
# plot 4:
check.summary <- cast(check.df, "MIxS + target.gene ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$target.gene
levels(check.summary$col) <- c("green", "yellow", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, 8000),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
barplot(
  as.matrix(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))]),
  col = levels(check.summary$col),
  legend.text = paste0(
    levels(check.summary$target.gene),
    " (", 
    round(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])/sum(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])) * 100),
    "%)"
  ),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
# plot 5:
check.df <- data.frame(
  nom.len = factor(ifelse(is.na(ena.out$nominal_length), "no", "yes"), levels = c("yes", "no")),
  year.created = gsub("-.*", "", ena.out$first_created),
  gfbio = ena.out$broker_name == "GFBIO",
  broker = ena.out$broker_name != "",
  study.accnos = ena.out$study_accession,
  MIxS = ena.out$environmental_package != "" & ena.out$checklist %in% mixs.checklists
)
check.summary <- cast(check.df, "nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, max(table(check.df$year.created))),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
abline(h = 8000)
barplot(
  as.matrix(check.summary[, -c(1, ncol(check.summary))]),
  col = as.character(check.summary$col),
  legend.text = paste0(check.summary$nom.len, " (", round(rowSums(check.summary[, -c(1, ncol(check.summary))])/sum(rowSums(check.summary[, -c(1, ncol(check.summary))])) * 100), "%)"),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
par(xpd = NA)
text(
  14.5, par("usr")[4],
  pos = 3,
  labels = "Nominal length",
  font = 2,
  cex = 1.2
)
par(xpd = F)
# plot 6:
check.summary <- cast(check.df, "MIxS + nom.len ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
check.summary$col <- check.summary$nom.len
levels(check.summary$col) <- c("green", "red")
plot(
  0, 0, 
  type = "n",
  axes = F,
  xlim = c(0.2, 13.2),
  ylim = c(0, 8000),
  lwd = 0.5
)
axis(
  2, 
  las = 2,
  cex.axis = 1,
  mgp = c(2.5, 0.5, 0),
  tcl = -0.3,
  lwd = 0.5
)
barplot(
  as.matrix(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))]),
  col = levels(check.summary$col),
  legend.text = paste0(
    levels(check.summary$nom.len), 
    " (", 
    round(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])/sum(rowSums(check.summary[as.logical(check.summary$MIxS), -c(1, 2, ncol(check.summary))])) * 100), 
    "%)"
  ),
  args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
  names.arg = colnames(check.summary)[-c(1, 2, ncol(check.summary))],
  add = T,
  axes = F,
  las = 2,
  mgp = c(2.5, 0.5, 0),
  cex.names = 1,
  lwd = 0.5
)
# labels
mtext("Submission year", 1, 0, outer = T)
mtext("Number of cases", 2, 0, outer = T)

dev.off()

### figure 2: envo terms ####

pdf("ENA_fig_envo.pdf", width = 7, height = 7.5)
# object size 700x750
par(mfcol = c(3, 2), mar = c(4, 4, 1, 0.5), oma = c(1.5, 2, 2, 0), ann = F, yaxs = "i", lwd = 0.5)

# plot 1-3:
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
check.summary <- list(
  cast(check.df, "biome ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T),
  cast(check.df, "material ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T),
  cast(check.df, "feature ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 1]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "blue", "red")
}
for(i in 1:length(check.summary)) {
  plot(
    0, 0, 
    type = "n",
    axes = F,
    xlim = c(0.2, 13.2),
    ylim = c(0, max(table(check.df$year.created))),
    lwd = 0.5
  )
  axis(
    2, 
    las = 2,
    cex.axis = 1,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3,
    lwd = 0.5
  )
  abline(h = 8000)
  barplot(
    as.matrix(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))]),
    col = as.character(check.summary[[i]]$col),
    legend.text = paste0(check.summary[[i]][, 1], " (", round(rowSums(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))])/sum(rowSums(check.summary[[i]][, -c(1, ncol(check.summary[[i]]))])) * 100), "%)"),
    args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
    names.arg = colnames(check.summary[[i]])[-c(1, ncol(check.summary[[i]]))],
    add = T,
    axes = F,
    las = 2,
    mgp = c(2.5, 0.5, 0),
    cex.names = 1,
    lwd = 0.5
  )
  if(i == 1) {
    mtext("All cases", 3, 1, outer = F)
  }
  par(xpd = NA)
  text(
    14.5, par("usr")[4],
    pos = 3,
    labels = names(check.summary)[i],
    font = 2,
    cex = 1.2
  )
  par(xpd = F)
}

# plot 4-6:
check.summary <- list(
  cast(check.df, "MIxS + biome ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T),
  cast(check.df, "MIxS + material ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T),
  cast(check.df, "MIxS + feature ~ year.created", value = "study.accnos", fun.aggregate = "length", fill = 0, add.missing = T)
)
names(check.summary) <- names(ena.in.envo)
for(i in 1:length(check.summary)) {
  check.summary[[i]]$col <- check.summary[[i]][, 2]
  levels(check.summary[[i]]$col) <- c("green", "yellow", "blue", "red")
}
for(i in 1:length(check.summary)) {
  plot(
    0, 0, 
    type = "n",
    axes = F,
    xlim = c(0.2, 13.2),
    ylim = c(0, 8000),
    lwd = 0.5
  )
  axis(
    2, 
    las = 2,
    cex.axis = 1,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3,
    lwd = 0.5
  )
  barplot(
    as.matrix(check.summary[[i]][as.logical(check.summary[[i]]$MIxS), -c(1, 2, ncol(check.summary[[i]]))]),
    col = levels(droplevels(check.summary[[i]][as.logical(check.summary[[i]]$MIxS), "col"])),
    legend.text = paste0(
      levels(droplevels(check.summary[[i]][as.logical(check.summary[[i]]$MIxS), 2])), 
      " (", 
      round(rowSums(check.summary[[i]][as.logical(check.summary[[i]]$MIxS), -c(1, 2, ncol(check.summary[[i]]))])/sum(rowSums(check.summary[[i]][as.logical(check.summary[[i]]$MIxS), -c(1, 2, ncol(check.summary[[i]]))])) * 100), 
      "%)"
    ),
    args.legend = list(x = "topleft", cex = 1.2, bty = "n", pt.cex = 1.5),
    names.arg = colnames(check.summary[[i]])[-c(1, 2, ncol(check.summary[[i]]))],
    add = T,
    axes = F,
    las = 2,
    mgp = c(2.5, 0.5, 0),
    cex.names = 1,
    lwd = 0.5
  )
  if(i == 1) {
    mtext("MIxS environmental package", 3, 1, outer = F)
  }
}
# labels
mtext("Submission year", 1, 0, outer = T)
mtext("Number of cases", 2, 0, outer = T)

dev.off()

### figure 3: venn diagram ####

venn(oligotyping.studies[, 2:4], ilcs = 1.3, sncs = 1)
