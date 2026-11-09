setwd("C:/Users/hassenru/Documents/Additional_projects/Data_mining_commentary/Update_2025/Summaries/")
require(tidyverse)
require(scales)
require(reshape)
require(patchwork)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
options(scipen=999)

### main figure 1: manually created in pptx ####


### main figure 2: general trends ####

x <- readRDS("summary_no_data.rds") %>% 
  ungroup()
plot_in <- x %>% 
  group_by(year_created) %>% 
  summarise(
    data_sets = sum(n_data), 
    data_volume = sum(mbp_data)/1000,
    ddbj_perc = sum(n_data[entry_db == "DRR"])/sum(n_data) * 100,
    ena_perc = sum(n_data[entry_db == "ERR"])/sum(n_data) * 100,
    ncbi_perc = sum(n_data[entry_db == "SRR"])/sum(n_data) * 100,
    mixs_perc = sum(n_data[mixs])/sum(n_data) * 100,
    gfbio_perc = sum(n_data[gfbio])/sum(n_data) * 100
  )
png("figure_general_trends.png", width = 7, height = 5, units = "in", res = 300)
par(mar = c(4, 3, 1, 8), xaxs = "i", xpd = NA)
bp <- barplot(
  rbind(
    plot_in$data_sets,
    rescale(plot_in$data_volume, from = c(0, max(plot_in$data_volume)), to = c(0, 650000))
  ),
  beside = T,
  axes = F,
  col = c("grey30", "grey70"),
  border = NA
)
axis(1, at = apply(bp, 2, mean), labels = plot_in$year_created, las = 2)
axis(4, at = seq(0, 650000, 100000), las = 3, col = "grey30", col.axis = "grey30", mgp = c(2, 0.4, 0.2), cex.axis = 0.8, tcl = -0.3)
mtext(text = "Number of runs", side = 4, line = 1.2, cex = 1, col = "grey30")
axis(4, at = rescale(seq(0, 30500, 5000), from = c(0, max(plot_in$data_volume)), to = c(0, 650000)), labels = seq(0, 30500, 5000), las = 3, col = "grey70", col.axis = "grey70", mgp = c(2, 3.4, 3.2), cex.axis = 0.8, tcl = -0.3)
mtext(text = "Data volume [GB]", side = 4, line = 4.2, cex = 1, col = "grey70")
for(i in 4:8) {
  lines(
    apply(bp, 2, mean),
    rescale(plot_in %>% pull(i), from = c(0, 100), to = c(0, 650000)), 
    col = c("darkred", "dodgerblue", "gold", "darkolivegreen", "magenta")[i - 3]
  )
  points(
    apply(bp, 2, mean), 
    rescale(plot_in %>% pull(i), from = c(0, 100), to = c(0, 650000)),
    bg = c("darkred", "dodgerblue", "gold", "darkolivegreen", "magenta")[i - 3],
    col = "white",
    pch = 21
  )
}
axis(2, at = rescale(seq(0, 100, 20), from = c(0, 100), to = c(0, 650000)), labels = seq(0, 100, 20), las = 3, col = "black", col.axis = "black", mgp = c(2, 0.5, 0.2), cex.axis = 0.8, tcl = -0.3)
mtext(text = "Percentage of runs", side = 2, line = 1.3, cex = 1, col = "black")
legend(
  "top",
  legend = c(
    paste0("DDBJ (", round(sum(x$n_data[x$entry_db == "DRR"])/sum(x$n_data) * 100),"%)"),
    paste0("ENA (", round(sum(x$n_data[x$entry_db == "ERR"])/sum(x$n_data) * 100),"%)"),
    paste0("NCBI (", round(sum(x$n_data[x$entry_db == "SRR"])/sum(x$n_data) * 100),"%)"),
    paste0("MIxS (", round(sum(x$n_data[x$mixs])/sum(x$n_data) * 100),"%)"),
    paste0("GFBio (", round(sum(x$n_data[x$gfbio])/sum(x$n_data) * 100, 1),"%)")
  ),
  col = c("darkred", "dodgerblue", "gold", "darkolivegreen", "magenta"),
  pch = 16,
  cex = 0.8,
  pt.cex = 1,
  bty = "n",
  ncol = 2
)
dev.off()


### main figure 3: geo location and date ####

x_files <- c(
  "summary_collection_date_available.rds",
  "summary_lat_lon_available.rds",
  "summary_z_coordinate.rds",
  "summary_country_available.rds"
)
names(x_files) <- gsub(".rds", "", gsub("summary_", "", x_files))
nice_name <- c("collection_date", "latitude/longitude", "vertical coordinates", "country")
names(nice_name) <- names(x_files)

x_col <-  c("grey", "darkred", "red", "orange", "dodgerblue", "blue", "navy")

x_bar <- lapply(
  x_files, 
  function(f){
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    levels(tmp$year_created) <- c("<=2014", "<=2014", "<=2014", "<=2014", "<=2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")
    levels(tmp[, 4])[4] <- "xml_syntax"
    levels(tmp[, 4])[7] <- "tsv_syntax"
    return(tmp)
  }
)

x_heat <- lapply(
  x_files, 
  function(f) {
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    levels(tmp[, 4])[4] <- "xml_syntax"
    levels(tmp[, 4])[7] <- "tsv_syntax"
    return(tmp)
  }
)

panel_id <- paste0("(", LETTERS[1:length(x_files)], ")")
names(panel_id) <- names(x_files)

png("figure_xyzt.png", height = 9.5, width = 11, unit = "in", res = 300)
layout(
  matrix(
    c(
      1, 2, 9, 10,
      3, 4, 11, 12,
      5, 6, 13, 14,
      7, 8, 15, 16
    ),
    nrow = 4,
    ncol = 4, 
    byrow = T
  )
)
par(ann = F, xpd = NA, lwd = 0.5, oma = c(3, 2, 3, 2))

for(i in names(x_files)) {
  par(mar = c(1, 0.5, 0, 2), xaxs = "r", yaxs = "r")
  tmp1 <- x_bar[[i]] %>% 
    filter(!mixs)  %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_bar[[i]] %>% 
    filter(mixs) %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  barplot(
    as.matrix(tmp1), 
    horiz = T,
    xlim = c(400000, 0), 
    names.arg = rep("", ncol(tmp1)),
    col = x_col,
    legend.text = paste0(rownames(tmp1), " (", round(rowSums(tmp1)/sum(tmp1) * 100), "%)"),
    args.legend = list(x = "topleft", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "collection_date_available") {
    mtext("without MIxS", side = 3, line = 0.5)
  }
  if(i == "country_available") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
    mtext("Number of runs", side = 1, line = 2.5, at = -50000, ad = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
  mtext(nice_name[i], side = 2, line = 0.5)
  text(
    par("usr")[1],
    par("usr")[4],
    pos = 2,
    labels = panel_id[i],
    cex = 1.2,
    font = 2
  )
  
  par(mar = c(1, 2, 0, 0.5))
  barplot(
    as.matrix(tmp2), 
    horiz = T, 
    xlim = c(0, 400000), 
    col = x_col,
    las = 1,
    legend.text = paste0(rownames(tmp2), " (", round(rowSums(tmp2)/sum(tmp2) * 100), "%)"),
    args.legend = list(x = "topright", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "collection_date_available") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
  if(i == "country_available") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
}

for(i in names(x_files)) {
  tmp1 <- x_heat[[i]] %>% 
    filter(!mixs)  %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_heat[[i]] %>% 
    filter(mixs) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp1))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 1, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
  title(ylab = "Percentage of runs", mgp = c(2, 0.5, 0))
  if(i == "country_available") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = colnames(tmp1), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp1)) {
    lines(
      (prop.table(t(tmp1), 1) * 100)[, j],
      lty = 1,
      col = x_col[j]
    )
    points(
      (prop.table(t(tmp1), 1) * 100)[, j],
      col = x_col[j],
      pch = 16
    )
  }
  if(i == "collection_date_available") {
    mtext("without MIxS", side = 3, line = 0.5)
  } 
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp2))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, labels = NA, lwd = 0.5)
  if(i == "country_available") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = colnames(tmp2), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp2)) {
    lines(
      (prop.table(t(tmp2), 1) * 100)[, j],
      lty = 1,
      col = x_col[j]
    )
    points(
      (prop.table(t(tmp2), 1) * 100)[, j],
      col = x_col[j],
      pch = 16
    )
  }
  if(i == "collection_date_available") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
}

dev.off()


### main figure 4: environment ####

x_files <- c(
  "summary_check_env_broad.rds",
  "summary_check_env_local.rds",
  "summary_check_env_medium.rds",
  "summary_env_tag_confidence.rds"
)
names(x_files) <- gsub(".rds", "", gsub("summary_", "", x_files))
nice_name <- c("broad-scale environment", "local environment", "environmental medium", "tag")
names(nice_name) <- names(x_files)

x_col <-  c("grey", "darkred", "orange", "dodgerblue", "blue")

x_bar <- lapply(
  x_files, 
  function(f){
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    levels(tmp$year_created) <- c("<=2014", "<=2014", "<=2014", "<=2014", "<=2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")
    return(tmp)
  }
)
levels(x_bar$env_tag_confidence$env_tag_confidence)[1] <- "no tag"

x_heat <- lapply(
  x_files, 
  function(f) {
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    return(tmp)
  }
)
levels(x_heat$env_tag_confidence$env_tag_confidence)[1] <- "no tag"

panel_id <- paste0("(", LETTERS[1:length(x_files)], ")")
names(panel_id) <- names(x_files)

png("figure_env.png", height = 9.5, width = 11, unit = "in", res = 300)
layout(
  matrix(
    c(
      1, 2, 9, 10,
      3, 4, 11, 12,
      5, 6, 13, 14,
      7, 8, 15, 16
    ),
    nrow = 4,
    ncol = 4, 
    byrow = T
  )
)
par(ann = F, xpd = NA, lwd = 0.5, oma = c(3, 2, 3, 2))

for(i in names(x_files)) {
  par(mar = c(1, 0.5, 0, 2), xaxs = "r", yaxs = "r")
  tmp1 <- x_bar[[i]] %>% 
    filter(!mixs)  %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_bar[[i]] %>% 
    filter(mixs) %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  barplot(
    as.matrix(tmp1), 
    horiz = T,
    xlim = c(400000, 0), 
    names.arg = rep("", ncol(tmp1)),
    col = x_col,
    legend.text = paste0(rownames(tmp1), " (", round(rowSums(tmp1)/sum(tmp1) * 100), "%)"),
    args.legend = list(x = "topleft", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "check_env_broad") {
    mtext("without MIxS", side = 3, line = 0.5)
  }
  if(i == "env_tag_confidence") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
    mtext("Number of runs", side = 1, line = 2.5, at = -50000, ad = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
  mtext(nice_name[i], side = 2, line = 0.5)
  text(
    par("usr")[1],
    par("usr")[4],
    pos = 2,
    labels = panel_id[i],
    cex = 1.2,
    font = 2
  )
  
  par(mar = c(1, 2, 0, 0.5))
  barplot(
    as.matrix(tmp2), 
    horiz = T, 
    xlim = c(0, 400000), 
    col = x_col,
    las = 1,
    legend.text = paste0(rownames(tmp2), " (", round(rowSums(tmp2)/sum(tmp2) * 100), "%)"),
    args.legend = list(x = "topright", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "check_env_broad") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
  if(i == "env_tag_confidence") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
}

for(i in names(x_files)) {
  tmp1 <- x_heat[[i]] %>% 
    filter(!mixs)  %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_heat[[i]] %>% 
    filter(mixs) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp1))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 1, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
  title(ylab = "Percentage of runs", mgp = c(2, 0.5, 0))
  if(i == "env_tag_confidence") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = colnames(tmp1), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp1)) {
    lines(
      (prop.table(t(tmp1), 1) * 100)[, j],
      lty = 1,
      col = x_col[j]
    )
    points(
      (prop.table(t(tmp1), 1) * 100)[, j],
      col = x_col[j],
      pch = 16
    )
  }
  if(i == "check_env_broad") {
    mtext("without MIxS", side = 3, line = 0.5)
  } 
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp2))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, labels = NA, lwd = 0.5)
  if(i == "env_tag_confidence") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = colnames(tmp2), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp2)) {
    lines(
      (prop.table(t(tmp2), 1) * 100)[, j],
      lty = 1,
      col = x_col[j]
    )
    points(
      (prop.table(t(tmp2), 1) * 100)[, j],
      col = x_col[j],
      pch = 16
    )
  }
  if(i == "check_env_broad") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
}

dev.off()


### main figure 5: marker gene, inconsistencies, publications ####

x_files <- c(
  "summary_target_gene_info.rds",
  "summary_metadata_inconsistencies.rds",
  "summary_publication_available.rds"
)
names(x_files) <- gsub(".rds", "", gsub("summary_", "", x_files))
nice_name <- c("amplicon", "inconsistencies", "scientific article")
names(nice_name) <- names(x_files)

x_col <- list(
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("darkred", "dodgerblue"),
  c("grey", "orange", "dodgerblue", "blue")
)
names(x_col) <- names(x_files)

x_bar <- lapply(
  x_files, 
  function(f){
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    levels(tmp$year_created) <- c("<=2014", "<=2014", "<=2014", "<=2014", "<=2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")
    return(tmp)
  }
)

x_heat <- lapply(
  x_files, 
  function(f) {
    tmp <- readRDS(f) %>% ungroup() %>% as.data.frame()
    return(tmp)
  }
)

panel_id <- paste0("(", LETTERS[1:length(x_files)], ")")
names(panel_id) <- names(x_files)

png("figure_gene.png", height = 7.5, width = 11, unit = "in", res = 300)
layout(
  matrix(
    c(
      1, 2, 7, 8,
      3, 4, 9, 10,
      5, 6, 11, 12
    ),
    nrow = 3,
    ncol = 4, 
    byrow = T
  )
)
par(ann = F, xpd = NA, lwd = 0.5, oma = c(3, 2, 3, 2))

for(i in names(x_files)) {
  par(mar = c(1, 0.5, 0, 2), xaxs = "r", yaxs = "r")
  tmp1 <- x_bar[[i]] %>% 
    filter(!mixs)  %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_bar[[i]] %>% 
    filter(mixs) %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  barplot(
    as.matrix(tmp1), 
    horiz = T,
    xlim = c(400000, 0), 
    names.arg = rep("", ncol(tmp1)),
    col = x_col[[i]],
    legend.text = paste0(rownames(tmp1), " (", round(rowSums(tmp1)/sum(tmp1) * 100), "%)"),
    args.legend = list(x = "topleft", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "target_gene_info") {
    mtext("without MIxS", side = 3, line = 0.5)
  }
  if(i == "publication_available") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
    mtext("Number of runs", side = 1, line = 2.5, at = -50000, ad = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
  mtext(nice_name[i], side = 2, line = 0.5)
  text(
    par("usr")[1],
    par("usr")[4],
    pos = 2,
    labels = panel_id[i],
    cex = 1.2,
    font = 2
  )
  
  par(mar = c(1, 2, 0, 0.5))
  barplot(
    as.matrix(tmp2), 
    horiz = T, 
    xlim = c(0, 400000), 
    col = x_col[[i]],
    las = 1,
    legend.text = paste0(rownames(tmp2), " (", round(rowSums(tmp2)/sum(tmp2) * 100), "%)"),
    args.legend = list(x = "topright", bty = "n", cex = 1, pt.cex = 1),
    axes = F
  )
  if(i == "target_gene_info") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
  if(i == "publication_available") {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
  } else {
    axis(1, cex.axis = 1, at = seq(0, 400000, 100000), labels = NA, lwd = 0.5)
  }
}

for(i in names(x_files)) {
  tmp1 <- x_heat[[i]] %>% 
    filter(!mixs)  %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x_heat[[i]] %>% 
    filter(mixs) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp1))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 1, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
  title(ylab = "Percentage of runs", mgp = c(2, 0.5, 0))
  if(i == "publication_available") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = colnames(tmp1), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp1)) {
    lines(
      (prop.table(t(tmp1), 1) * 100)[, j],
      lty = 1,
      col = x_col[[i]][j]
    )
    points(
      (prop.table(t(tmp1), 1) * 100)[, j],
      col = x_col[[i]][j],
      pch = 16
    )
  }
  if(i == "target_gene_info") {
    mtext("without MIxS", side = 3, line = 0.5)
  } 
  
  par(mar = c(1, 2, 0, 0))
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp2))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, labels = NA, lwd = 0.5)
  if(i == "publication_available") {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = colnames(tmp2), lwd = 0.5, las = 2)
  } else {
    axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = NA, lwd = 0.5)
  }
  for(j in 1:nrow(tmp2)) {
    lines(
      (prop.table(t(tmp2), 1) * 100)[, j],
      lty = 1,
      col = x_col[[i]][j]
    )
    points(
      (prop.table(t(tmp2), 1) * 100)[, j],
      col = x_col[[i]][j],
      pch = 16
    )
  }
  if(i == "target_gene_info") {
    mtext("with MIxS", side = 3, line = 0.5)
  }
}

dev.off()


### main figure 6: overall quality ####

x_files <- "summary_overall_quality.rds"
x_col <- c("darkred", "orange", "dodgerblue")

x_bar <- readRDS(x_files) %>% ungroup() %>% as.data.frame()
levels(x_bar$year_created) <- c("<=2014", "<=2014", "<=2014", "<=2014", "<=2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")

x_heat <- readRDS(x_files) %>% ungroup() %>% as.data.frame()

png("figure_overall_quality.png", height = 4, width = 11, unit = "in", res = 300)
layout(
  matrix(
    c(
      1, 2, 3, 4
    ),
    nrow = 1,
    ncol = 4, 
    byrow = T
  )
)
par(ann = F, xpd = NA, lwd = 0.5, oma = c(3, 2, 3, 2))

par(mar = c(1, 0.5, 0, 2), xaxs = "r", yaxs = "r")
tmp1 <- x_bar %>% 
  filter(!mixs)  %>% 
  mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
  cast("overall_quality ~ year_created", value = "n_data", fun.aggregate = sum) %>% 
  as.data.frame()
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[, -1]
tmp2 <- x_bar %>% 
  filter(mixs) %>% 
  mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
  cast("overall_quality ~ year_created", value = "n_data", fun.aggregate = sum) %>% 
  as.data.frame()
rownames(tmp2) <- tmp2[, 1]
tmp2 <- tmp2[, -1]

barplot(
  as.matrix(tmp1), 
  horiz = T,
  xlim = c(400000, 0), 
  names.arg = rep("", ncol(tmp1)),
  col = x_col,
  legend.text = paste0(rownames(tmp1), " (", round(rowSums(tmp1)/sum(tmp1) * 100), "%)"),
  args.legend = list(x = "topleft", bty = "n", cex = 1, pt.cex = 1),
  axes = F
)
axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
mtext("Number of runs", side = 1, line = 2.5, at = -50000, ad = 0.5)
mtext("without MIxS", side = 3, line = 0.5)

par(mar = c(1, 2, 0, 0.5))
barplot(
  as.matrix(tmp2), 
  horiz = T, 
  xlim = c(0, 400000), 
  col = x_col,
  las = 1,
  legend.text = paste0(rownames(tmp2), " (", round(rowSums(tmp2)/sum(tmp2) * 100), "%)"),
  args.legend = list(x = "topright", bty = "n", cex = 1, pt.cex = 1),
  axes = F
)
axis(1, cex.axis = 1, at = seq(0, 400000, 100000), lwd = 0.5)
mtext("with MIxS", side = 3, line = 0.5)
 
tmp1 <- x_heat %>% 
  filter(!mixs)  %>% 
  cast("overall_quality ~ year_created", value = "n_data", fun.aggregate = sum) %>% 
  as.data.frame()
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[, -1]
tmp2 <- x_heat %>% 
  filter(mixs) %>% 
  cast("overall_quality ~ year_created", value = "n_data", fun.aggregate = sum) %>% 
  as.data.frame()
rownames(tmp2) <- tmp2[, 1]
tmp2 <- tmp2[, -1]

par(mar = c(1, 2, 0, 0))
plot(
  0, 0, 
  type = "n",
  axes = F,
  ylim = c(0, 100),
  xlim = c(1, ncol(tmp1))
)
box("plot")
axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 1, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
title(ylab = "Percentage of runs", mgp = c(2, 0.5, 0))
axis(1, cex.axis = 1, at = 1:ncol(tmp1), labels = colnames(tmp1), lwd = 0.5, las = 2)
for(j in 1:nrow(tmp1)) {
  lines(
    (prop.table(t(tmp1), 1) * 100)[, j],
    lty = 1,
    col = x_col[j]
  )
  points(
    (prop.table(t(tmp1), 1) * 100)[, j],
    col = x_col[j],
    pch = 16
  )
}
mtext("without MIxS", side = 3, line = 0.5)
 
par(mar = c(1, 2, 0, 0))
plot(
  0, 0, 
  type = "n",
  axes = F,
  ylim = c(0, 100),
  xlim = c(1, ncol(tmp2))
)
box("plot")
axis(2, at = seq(0, 100, 20), las = 2, labels = NA, lwd = 0.5)
axis(1, cex.axis = 1, at = 1:ncol(tmp2), labels = colnames(tmp2), lwd = 0.5, las = 2)
for(j in 1:nrow(tmp2)) {
  lines(
    (prop.table(t(tmp2), 1) * 100)[, j],
    lty = 1,
    col = x_col[j]
  )
  points(
    (prop.table(t(tmp2), 1) * 100)[, j],
    col = x_col[j],
    pch = 16
  )
}
mtext("with MIxS", side = 3, line = 0.5)

dev.off()


### main figure 7: case studies ####

# Molari
molari.usecase <- read.table("molari_sets.txt", h = T, sep = "\t")
p1 <- upset(
  molari.usecase,
  intersect = colnames(molari.usecase), 
  width_ratio = 0.4,
  name = "Molari et al. (2023)",
  base_annotations = list(
    'Intersection size' = (
      intersection_size() +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() +
        theme(
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    )
  ),
  set_sizes = (
    upset_set_size() +
      theme_classic() +
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  ),
  stripes='white'
) 

# Riedinger
riedinger.usecase <- read.table("riedinger_usecase.txt", h = T, sep = "\t") %>% 
  mutate(
    primer_clipped = as.numeric(primer_present == "removed"),
    primer_included = as.numeric(primer_present == "included"),
    merged = as.numeric(merged_status == "pe_merged"),
    unmerged = as.numeric(merged_status == "pe_separate")
  )
upset_in <- riedinger.usecase %>% select(primer_clipped, primer_included, merged, unmerged)
p2 <- upset(
  upset_in, 
  intersect = colnames(upset_in), 
  width_ratio = 0.4,
  name = "Riedinger et al. (2025)",
  base_annotations = list(
    'Intersection size' = (
      intersection_size() +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() +
        theme(
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    )
  ),
  set_sizes = (
    upset_set_size() +
      theme_classic() +
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  ),
  stripes='white'
)

# Svae plot
ggsave("case_studies.png", wrap_elements(p1) + wrap_elements(p2), width = 9, height = 5, units = "in", dpi = 300)



### SI figure 1: frequent taxids ####

x <- read.table("stats_taxid.txt", h = T, sep = "\t")
plot(x$n_data, x$mbp_data)
topN <- x$scientific_name[order(x$n_data, decreasing = T)[1:25]]
plot_in <- x %>% 
  mutate(
    tax_group = ifelse(
      scientific_name %in% topN,
      scientific_name,
      ifelse(
        grepl("metagenome", scientific_name),
        "other metagenomes",
        "other unclassified sequences"
      )
    )
  ) %>% 
  group_by(tax_group) %>% 
  summarise(
    data_sets = sum(n_data),
    data_volume = sum(mbp_data)/1000
  ) %>% 
  arrange(grepl("other", tax_group), desc(data_sets))
png("figure_taxids.png", width = 7, height = 5, units = "in", res = 300)
par(mar = c(8, 4, 1, 4), xpd = NA, xaxs = "i")
bp <- barplot(
  rbind(
    plot_in$data_sets,
    rescale(plot_in$data_volume, from = c(0, max(plot_in$data_volume)), to = c(0, 750000))
  ),
  beside = T,
  axes = F,
  col = c("grey30", "grey70"),
  border = NA
)
axis(1, at = apply(bp, 2, mean), labels = NA)
text(
  apply(bp, 2, mean),
  par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.04,
  labels = plot_in$tax_group,
  cex = 0.8,
  srt = 45,
  adj = 1
)
axis(2, at = seq(0, 750000, 100000), las = 2, col = "grey30", col.axis = "grey30", cex.axis = 0.8, mgp = c(2, 0.6, 0.2), tcl = -0.3)
mtext(text = "Number of runs", side = 2, line = 3, cex = 1, col = "grey30")
axis(4, at = rescale(seq(0, 25000, 5000), from = c(0, max(plot_in$data_volume)), to = c(0, 750000)), labels = seq(0, 25000, 5000), las = 2, col = "grey70", col.axis = "grey70", cex.axis = 0.8, mgp = c(2, 0.6, 0.2), tcl = -0.3)
mtext(text = "Data volume [GB]", side = 4, line = 2.8, cex = 1, col = "grey70")
dev.off()


### preliminary barplots ####

x_files <- c(
  "summary_collection_date_available.rds",
  "summary_lat_lon_available.rds",
  "summary_z_coordinate.rds",
  "summary_country_available.rds",
  "summary_check_env_broad.rds",
  "summary_check_env_local.rds",
  "summary_check_env_medium.rds",
  "summary_env_tag_confidence.rds",
  "summary_target_gene_info.rds",
  "summary_metadata_inconsistencies.rds",
  "summary_publication_available.rds",
  "summary_overall_quality.rds"
)
names(x_files) <- gsub(".rds", "", gsub("summary_", "", x_files))

x_col <- list(
  c("grey", "darkred", "red", "orange", "dodgerblue", "blue", "navy"),
  c("grey", "darkred", "red", "orange", "dodgerblue", "blue", "navy"),
  c("grey", "darkred", "red", "orange", "dodgerblue", "blue", "navy"),
  c("grey", "darkred", "red", "orange", "dodgerblue", "blue", "navy"),
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("grey", "darkred", "orange", "dodgerblue", "blue"),
  c("darkred", "dodgerblue"),
  c("grey", "orange", "dodgerblue", "blue"),
  c("darkred", "orange", "dodgerblue")
)
names(x_col) <- names(x_files)

x <- lapply(
  x_files, 
  function(f){
    tmp <- readRDS(f) %>% ungroup()
    levels(tmp$year_created) <- c("<=2014", "<=2014", "<=2014", "<=2014", "<=2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")
    return(tmp)
  }
)

pdf("figure_trends_barplot.pdf", height = 6, width = 7, onefile = T)
for(i in names(x_files)) {
  par(mfcol = c(1, 2), mar = c(3, 3, 3, 2), lwd = 0.5, oma = c(0, 0, 1, 0))
  tmp1 <- x[[i]] %>% 
    filter(!mixs)  %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x[[i]] %>% 
    filter(mixs) %>% 
    mutate(year_created = ordered(year_created, levels = rev(levels(year_created)))) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  barplot(
    as.matrix(tmp1), 
    horiz = T,
    xlim = c(400000, 0), 
    names.arg = rep("", ncol(tmp1)),
    col = x_col[[i]],
    legend.text = paste0(
      rownames(tmp1),
      " (",
      round(rowSums(tmp1)/sum(tmp1) * 100),
      "%)"
    ),
    args.legend = list(x = "topleft", bty = "n", cex = 0.8, pt.cex = 1),
    cex.axis = 0.8
  )
  mtext("no MIxS", side = 3, line = 2)
  par(mar = c(3, 2, 3, 3))
  barplot(
    as.matrix(tmp2), 
    horiz = T, 
    xlim = c(0, 400000), 
    col = x_col[[i]],
    las = 1,
    legend.text = paste0(
      rownames(tmp2),
      " (",
      round(rowSums(tmp2)/sum(tmp2) * 100),
      "%)"
    ),
    args.legend = list(x = "topright", bty = "n", cex = 0.8, pt.cex = 1),
    cex.axis = 0.8
  )
  mtext("with MIxS", side = 3, line = 2)
  mtext(i, outer = T, side = 3, line = 0)
}
dev.off()


### preliminary heatmaps ####

x <- lapply(x_files, function(f) readRDS(f) %>% ungroup())

pdf("figure_trends_heatmap.pdf", height = 4.5, width = 7, onefile = T)
for(i in names(x)) {
  tmp1 <- x[[i]] %>% 
    filter(!mixs)  %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp1) <- tmp1[, 1]
  tmp1 <- tmp1[, -1]
  tmp2 <- x[[i]] %>% 
    filter(mixs) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp2) <- tmp2[, 1]
  tmp2 <- tmp2[, -1]
  
  layout(
    matrix(1:4, nrow = 2, ncol = 2, byrow = T),
    widths = c(5, 1),
    heights = c(1, 1)
  )
  par(ann = F, mar = c(0.2, 9, 3, 0), xaxs = "i", yaxs = "i", xpd = NA)
  image(
    prop.table(t(tmp1), 1) * 100,
    col = colorRampPalette(c("grey95", "orange", "red", "darkred"))(50),
    y = 1:nrow(tmp1),
    x = 1:ncol(tmp1),
    zlim = c(0, 100),
    axes = F
  )
  box("plot")
  axis(2, at = 1:nrow(tmp1), labels = rownames(tmp1), las = 2, cex.axis = 0.8)
  mtext("no MIxS", side = 2, line = 8)
  mtext(i, side = 3, line = 1)
  
  par(mar = c(0.2, 1, 3, 4))
  image(
    matrix(
      seq(0, 100, length.out = 50),
      ncol = 50,
      nrow = 1
    ),
    y = seq(0, 100, length.out = 50),
    x = 1,
    axes = F,
    col = colorRampPalette(c("grey95", "orange", "red", "darkred"))(50)
  )
  box("plot")
  axis(
    4, 
    at = seq(0, 100, 20), 
    mgp = c(1, 0.5, 0),
    tcl = -0.3,
    cex.axis = 0.8,
    las = 2
  )
  
  par(mar = c(3, 9, 0.2, 0))
  image(
    prop.table(t(tmp2), 1) * 100,
    col = colorRampPalette(c("grey95", "dodgerblue", "blue", "navy"))(50),
    y = 1:nrow(tmp2),
    x = 1:ncol(tmp2),
    zlim = c(0, 100),
    axes = F
  )
  box("plot")
  axis(2, at = 1:nrow(tmp2), labels = rownames(tmp2), las = 2, cex.axis = 0.8)
  mtext("with MIxS", side = 2, line = 8)
  axis(1, at = 1:ncol(tmp2), labels = colnames(tmp2), las = 2)
  
  par(mar = c(3, 1, 0.2, 4))
  image(
    matrix(
      seq(0, 100, length.out = 50),
      ncol = 50,
      nrow = 1
    ),
    y = seq(0, 100, length.out = 50),
    x = 1,
    axes = F,
    col = colorRampPalette(c("grey95", "dodgerblue", "blue", "navy"))(50)
  )
  box("plot")
  axis(
    4, 
    at = seq(0, 100, 20), 
    mgp = c(1, 0.5, 0),
    tcl = -0.3,
    cex.axis = 0.8,
    las = 2
  )
}
dev.off()


### preliminary entry DB ####

x <- lapply(x_files, function(f) readRDS(f) %>% ungroup())

pdf("figure_entry_db_barplot.pdf", height = 5, width = 7, onefile = T)
for(i in names(x_files)) {
  par(mfcol = c(1, 3), mar = c(3, 2, 3, 1), lwd = 0.5, oma = c(0, 0, 2, 0), lwd = 0.5)
  for(j in unique(x[[i]]$entry_db)) {
    tmp <- x[[i]] %>% 
      filter(entry_db == j)  %>% 
      cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
      as.data.frame()
    rownames(tmp) <- tmp[, 1]
    tmp <- tmp[, -1]
    barplot(
      as.matrix(tmp), 
      horiz = F,
      col = x_col[[i]],
      legend.text = paste0(
        rownames(tmp),
        " (",
        round(rowSums(tmp)/sum(tmp) * 100),
        "%)"
      ),
      args.legend = list(x = "topleft", bty = "n", cex = 0.8, pt.cex = 1),
      cex.axis = 0.8,
      las = 2,
      axes = F
    )
    axis(2, las = 3, lwd = 0.5, cex.axis = 0.8, tcl = -0.3, mgp = c(2, 0.5, 0))
    mtext(j, side = 3, line = 1)
  }
  mtext(i, outer = T, side = 3, line = 0)
}
dev.off()

pdf("figure_entry_db_lineplot.pdf", height = 5, width = 7, onefile = T)
for(i in names(x_files)) {
  par(mfrow = c(1, 3), oma = c(3, 3, 3, 0), mar = c(1, 1, 3, 1))
  for(j in unique(x[[i]]$entry_db)) {
    tmp <- x[[i]] %>% 
      filter(entry_db == j)  %>% 
      cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
      as.data.frame()
    rownames(tmp) <- tmp[, 1]
    tmp <- tmp[, -1]
    plot(
      0, 0, 
      type = "n",
      axes = F,
      ylim = c(0, 100),
      xlim = c(1, ncol(tmp))
    )
    box("plot")
    if(j == unique(x[[i]]$entry_db)[1]) {
      axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 0.8, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
      mtext(text = "Percentage of runs", outer = T, line = 1, side = 2)
    } else {
      axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 0.8, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3, labels = NA)
    }
    for(k in 1:nrow(tmp)) {
      lines(
        (prop.table(t(tmp), 1) * 100)[, k],
        lty = 1,
        col = x_col[[i]][k]
      )
      points(
        (prop.table(t(tmp), 1) * 100)[, k],
        col = x_col[[i]][k],
        pch = 16
      )
    }
    mtext(text = j, side = 3, line = 1)
    axis(1, cex.axis = 0.8, at = 1:ncol(tmp), labels = colnames(tmp), lwd = 0.5, las = 2)
    mtext(text = i, side = 3, line = 1, outer = T)
  }
}
dev.off()

col_heat <- list(
  c("grey95", "orange", "red", "darkred"),
  c("grey95", "dodgerblue", "blue", "navy"),
  c("grey95", "greenyellow", "chartreuse3", "darkgreen")
)
names(col_heat) <- unique(x[[1]]$entry_db)
pdf("figure_entry_db_heatmap.pdf", height = 4.5, width = 7, onefile = T)
for(i in names(x_files)) {
  layout(
    matrix(1:6, nrow = 3, ncol = 2, byrow = T),
    widths = c(5, 1),
    heights = c(1, 1, 1)
  )
  par(oma = c(3, 0, 2, 0))
  for(j in unique(x[[i]]$entry_db)) {
    tmp <- x[[i]] %>% 
      filter(entry_db == j)  %>% 
      cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
      as.data.frame()
    rownames(tmp) <- tmp[, 1]
    tmp <- tmp[, -1]
    par(ann = F, mar = c(0.2, 9, 0.2, 0), xaxs = "i", yaxs = "i", xpd = NA)
    image(
      prop.table(t(tmp), 1) * 100,
      col = colorRampPalette(col_heat[[j]])(50),
      y = 1:nrow(tmp),
      x = 1:ncol(tmp),
      zlim = c(0, 100),
      axes = F
    )
    box("plot")
    axis(2, at = 1:nrow(tmp), labels = rownames(tmp), las = 2, cex.axis = 0.8)
    mtext(j, side = 2, line = 7.5)
    if(j == "SRR") {
      axis(1, at = 1:ncol(tmp), labels = colnames(tmp), las = 2)
    }
    par(mar = c(0.2, 1, 0.2, 4))
    image(
      matrix(
        seq(0, 100, length.out = 50),
        ncol = 50,
        nrow = 1
      ),
      y = seq(0, 100, length.out = 50),
      x = 1,
      axes = F,
      col = colorRampPalette(col_heat[[j]])(50)
    )
    box("plot")
    axis(
      4, 
      at = seq(0, 100, 20), 
      mgp = c(1, 0.5, 0),
      tcl = -0.3,
      cex.axis = 0.8,
      las = 2
    )
  }
  mtext(i, outer = T, side = 3, line = 0.5)
}
dev.off()


### preliminary GFBio ####

x <- lapply(x_files, function(f) readRDS(f) %>% ungroup())

pdf("figure_gfbio_barplot.pdf", height = 6, width = 7, onefile = T)
for(i in names(x)) {
  tmp <- x[[i]] %>% 
    filter(gfbio)  %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp) <- tmp[, 1]
  tmp <- tmp[, -1]
  tmp <- tmp[, colSums(tmp) > 1]
  barplot(
    as.matrix(tmp), 
    horiz = F,
    col = x_col[[i]],
    legend.text = paste0(
      rownames(tmp),
      " (",
      round(rowSums(tmp)/sum(tmp) * 100),
      "%)"
    ),
    args.legend = list(x = "topleft", bty = "n", cex = 0.8, pt.cex = 1),
    cex.axis = 0.8,
    las = 2,
    axes = F
  )
  axis(2, las = 3, lwd = 0.5, cex.axis = 0.8, tcl = -0.3, mgp = c(2, 0.5, 0))
  mtext(i, side = 3, line = 1)
}
dev.off()

pdf("figure_gfbio_lineplot.pdf", height = 5, width = 5, onefile = T)
for(i in names(x_files)) {
  par(oma = c(3, 3, 3, 0), mar = c(1, 1, 3, 1))
  tmp <- x[[i]] %>% 
    filter(gfbio) %>% 
    cast(paste0(i, " ~ year_created"), value = "n_data", fun.aggregate = sum) %>% 
    as.data.frame()
  rownames(tmp) <- tmp[, 1]
  tmp <- tmp[, -1]
  plot(
    0, 0, 
    type = "n",
    axes = F,
    ylim = c(0, 100),
    xlim = c(1, ncol(tmp))
  )
  box("plot")
  axis(2, at = seq(0, 100, 20), las = 2, cex.axis = 0.8, lwd = 0.5, mgp = c(2, 0.5, 0), tcl = -0.3)
  mtext(text = "Percentage of runs", outer = T, line = 1, side = 2)
  for(k in 1:nrow(tmp)) {
    lines(
      (prop.table(t(tmp), 1) * 100)[, k],
      lty = 1,
      col = x_col[[i]][k]
    )
    points(
      (prop.table(t(tmp), 1) * 100)[, k],
      col = x_col[[i]][k],
      pch = 16
    )
  }
  axis(1, cex.axis = 0.8, at = 1:ncol(tmp), labels = colnames(tmp), lwd = 0.5, las = 2)
  mtext(text = i, side = 3, line = 1, outer = T)
}
dev.off()


# only 124 studies brokered by GFBio
# this is likely correct
### preliminary world map ####


require(maps)
xy <- ena.assess %>% select(lon, lat, overall_quality) %>% as.data.frame()
xy$col <- xy$overall_quality
levels(xy$col) <- c("darkred", "orange", "dodgerblue")
xy <- xy[order(xy$overall_quality), ]

png("Summaries/world_map_sampling_locations.png", width = 15, height = 10, units = "in", res = 600)
maps::map(
  "world", 
  resolution = 2,
  fill = T, 
  col = "grey90",
  border = NA,
  mar = c(1, 1, 1, 1)
)
box("plot")
points(
  xy$lon,
  xy$lat,
  pch = 16,
  cex = ifelse(xy$overall_quality == "machine_readable", 0.8, 0.3),
  col = as.character(xy$col)
)
dev.off()

# no color
png("Summaries/world_map_sampling_locations_black.png", width = 15, height = 10, units = "in", res = 600)
maps::map(
  "world", 
  resolution = 2,
  fill = T, 
  col = "grey90",
  border = NA,
  mar = c(1, 1, 1, 1)
)
box("plot")
points(
  xy$lon,
  xy$lat,
  pch = 16,
  cex = 0.5,
  col = "black"
)
dev.off()
