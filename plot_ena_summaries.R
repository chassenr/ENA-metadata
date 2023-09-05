# plot ENA metadata status

plot_ena_summaries <- function(
    dat_all,
    dat_mixs,
    col,
    ymax_all = 180000,
    ymax_mixs = 80000
) {
  
  # format x axis
  bp <- barplot(
    as.matrix(
      data.frame(dat_all[, -1], check.names = F)
    ),
    plot = F
  )
  
  # plot 1: complete data set
  plot(
    0, 0, 
    type = "n",
    axes = F,
    xlim = c(min(bp), max(bp)),
    ylim = c(0, ymax_all),
    lwd = 0.5
  )
  axis(
    2, 
    at = seq(0, ymax_all, 20000),
    las = 2,
    cex.axis = 1.3,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3,
    lwd = 0.5
  )
  # abline(h = ymax_mixs)
  barplot(
    as.matrix(
      data.frame(dat_all[, -1], check.names = F)
    ),
    col = col,
    legend.text = dat_all[, 1],
    args.legend = list(
      x = "topleft", 
      bty = "n",
      # legend = rev(paste0(dat_all[, 1], " (", round(rowSums(dat_all[, -1])/sum(dat_all[, -1]) * 100), "%)")),
      cex = 1.3
    ),
    names.arg = rep("", length(bp)),
    add = T,
    axes = F,
    lwd = 0.5
  )
  par(xpd = NA)
  tmp_prop <- prop.table(
    as.matrix(
      data.frame(dat_all[, -1], check.names = F)
    ),
    2
  ) * 100
  tmp_prop <- tmp_prop[nrow(tmp_prop):1, ]
  text(
    rep(bp, each = nrow(dat_all)),
    rep(line2user(seq(0, 2.5, length.out = nrow(dat_all)), side = 1), length(bp)),
    # pos = 1,
    labels = ifelse(is.na(tmp_prop), "", round(tmp_prop)),
    cex = 0.8
  )
  text(
    rep(line2user(0, side = 2), nrow(dat_all)),
    line2user(seq(0, 2.5, length.out = nrow(dat_all)), side = 1),
    # pos = 1,
    labels = rev(dat_all[, 1]),
    adj = 1,
    cex = 0.8
  )
  par(xpd = F)
  axis(1, at = bp, col = NA, col.ticks = NA, labels = colnames(dat_all)[-1], las = 2, mgp = c(3, 3, 0), hadj = 1, cex.axis = 1.3)
  
  # plot 2: MIxS checklists used
  plot(
    0, 0, 
    type = "n",
    axes = F,
    xlim = c(min(bp), max(bp)),
    ylim = c(0, ymax_mixs),
    lwd = 0.5
  )
  axis(
    2, 
    at = seq(0, ymax_mixs, 10000),
    las = 2,
    cex.axis = 1.3,
    mgp = c(2.5, 0.5, 0),
    tcl = -0.3,
    lwd = 0.5
  )
  barplot(
    as.matrix(
      data.frame(dat_mixs[, -1], check.names = F)
    ),
    ylim = c(0, ymax_mixs),
    names.arg = rep("", length(bp)),
    col = col,
    add = T,
    axes = F,
    lwd = 0.5
    # legend.text = dat_mixs[, 1],
    # args.legend = list(
    #   x = "topleft",
    #   # legend = rev(paste0(dat_mixs[, 1], " (", round(rowSums(dat_mixs[, -1])/sum(dat_mixs[, -1]) * 100), "%)"))
    #   bty = "n",
    #   cex = 1.3
    # )
  )
  par(xpd = NA)
  tmp_prop <- prop.table(
    as.matrix(
      data.frame(dat_mixs[, -1], check.names = F)
    ),
    2
  ) * 100
  tmp_prop <- tmp_prop[nrow(tmp_prop):1, ]
  text(
    rep(bp, each = nrow(dat_mixs)),
    rep(line2user(seq(0, 2.5, length.out = nrow(dat_mixs)), side = 1), length(bp)),
    # pos = 1,
    labels = ifelse(is.na(tmp_prop), "", round(tmp_prop)),
    cex = 0.8
  )
  # text(
  #   rep(line2user(0, side = 2), nrow(dat_mixs)),
  #   line2user(seq(0, 2.5, length.out = nrow(dat_mixs)), side = 1),
  #   # pos = 1,
  #   labels = rev(dat_mixs[, 1]),
  #   adj = 1,
  #   cex = 0.8
  # )
  par(xpd = F)
  axis(1, at = bp, col = NA, col.ticks = NA, labels = colnames(dat_mixs)[-1], las = 2, mgp = c(3, 3, 0), hadj = 1, cex.axis = 1.3)
  
}
