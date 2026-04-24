# Import packages
library(Biostrings)
library(brms)
library(coda)
library(ggtree)
library(MASS)
library(tidyverse)
library(ggridges)
library(ggcorrplot)
library(gridExtra)
library(gtools)
library(phytools)
library(stringr)
library(vioplot)

# Functions
######
# show posterior predictions as a function of depth,
# plot observed data on top to compare
depth_ppd <- function(model, data, file_out, depth_metric, main, norm, mean = FALSE) {
  pdf(file_out, width = 5, height = 5)
  par(xpd = FALSE)
  trans_light_blue <- rgb(135, 206, 250, alpha = 20, 
                    maxColorValue = 255)
  # get posterior prediction
  ppd <- posterior_predict(model)
  # if you want to show posterior means instead
  if (mean == TRUE) {
    ppd <- posterior_epred(model,
                           re_formula = NA)
  }
  maximum <- 1000
  # whether you want to plot out against normalised depth
  # or the raw depth
  xlab <- "Normalised depth"
  if (norm == FALSE) {
    xlab = "Depth"
  }
  # just to prepare the plotting area
  plot(ppd[1,] ~ data[, depth_metric],
       ylim = c(0, maximum), pch = "",
       xlab = xlab,
       ylab = "# of integrations", main = main)
  # plot PPD
  for (col in 1:dim(ppd)[2]) {
    points(ppd[, col] ~ rep(data[, depth_metric][col], dim(ppd)[1]),
           col = trans_light_blue, pch = 20)
  }
  # plot observed data
  points(data$int ~ data[, depth_metric], pch = 20, col = "steelblue")
  dev.off()
}

# print out the 95% HPDI for the difference between the 
# expected number of integrations
# at two different depths
depth_effect_HPDI <- function(model, depths, original_vector) {
  # get mean and SD of original vector for scaling
  scale_mean <- mean(original_vector)
  scale_sd <- sd(original_vector)
  # extract draws
  draws <- as.matrix(model)
  # get the coefficient for depth
  depth_coefs <- draws[,2]
  # get the intercept
  intercepts <- draws[,1]
  # convert the provided raw depth values to the cube root
  transf <- depths[1] ^ (1/3)
  print(transf ^ 3)
  transf <- (transf - scale_mean)/scale_sd
  # get expected counts at first depth value
  epred1 <- exp(intercepts + transf * depth_coefs)
  transf <- depths[2] ^ (1/3)
  print(transf ^ 3)
  transf <- (transf - scale_mean)/scale_sd
  # get expected counts at second depth value
  epred2 <- exp(intercepts + transf * depth_coefs)
  # get difference 
  diffs <- epred2 - epred1
  print("median of the posterior:")
  print(median(diffs))
  print("95% HPDI:")
  print(HPDinterval(as.mcmc(diffs), prob = 0.95))
  print("% above 0:")
  print(mean(diffs > 0) * 100)
}

plot_order_corrs <- function(depth_ordered, dataset, wasps, parameter = "depth_cube_root_cpm", filter = c(),
                             highlight = c()) {
  depth_ordered <- as.character(depth_ordered)
  if (length(filter) > 0) {
    depth_ordered_copy <- depth_ordered
    depth_ordered <- depth_ordered[depth_ordered %in% filter]
  }
  depths <- matrix(0, nrow = length(depth_ordered),
                   ncol = length(wasps))
  circles_clean <- rep(0, length(depth_ordered))
  for (seg in 1:length(depth_ordered)) {
    curr_seg <- str_split(depth_ordered[seg], "_")[[1]]
    if (length(curr_seg) > 2) {
      curr_seg <- paste(curr_seg[2], curr_seg[3], sep = "_")
    } else {
      curr_seg <- curr_seg[2]
    }
    circles_clean[seg] <- curr_seg
  }
  wasps <- rev(c("C. sesamiae kitale", "C. sesamiae mombasa",
             "C. typhae", "C. flavipes", "C. icipe"))
  for (wasp in 1:length(wasps)) {
    curr_dataset <- dataset[dataset$wasp == wasps[wasp],]
    agg <- aggregate(curr_dataset[, parameter] ~ curr_dataset[, "Segment"],
                     FUN = mean)
    rownames(agg) <- agg[,1]
    depths[, wasp] <- agg[depth_ordered, 2]
  }
  par(mfrow = c(1,1))
  cor_matrix <- cor(depths, method = "spearman")
  wasps_temp <- as.character(wasps)
  wasps_temp <- lapply(wasps_temp, FUN = function(x) {y <- strsplit(x, " ")[[1]]; y[length(y)]})
  wasps_temp <- rev(c("C. s. Kitale", "C. s. Mombasa", "C. typhae", "C. flavipes", "C. icipe"))
  colnames(depths) <- wasps_temp
  rownames(cor_matrix) <- wasps_temp
  colnames(cor_matrix) <- wasps_temp
  all_p <- c()
  seen <- c()
  for (w1 in wasps_temp) {
    for (w2 in wasps_temp) {
      if (w1 != w2) {
        possible_strings <- c(paste(w1, w2), paste(w2, w1))
        if (!(possible_strings[1] %in% seen) & !(possible_strings[2] %in% seen)) {
          seen <- c(seen, possible_strings[1])
          p <- cor.test(depths[, w1], depths[, w2], method = "spearman")$p.value
          all_p <- c(all_p, p)
        }
      }
    }
  }
  all_p <- signif(p.adjust(all_p, method = "holm"), 1)
  p_mat <- matrix(1000, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  colnames(p_mat) <- wasps_temp
  rownames(p_mat) <- wasps_temp
  for (w1 in wasps_temp) {
    for (w2 in wasps_temp) {
      if (w1 == w2) {
        p_mat[w1, w2] <- 0
      } else {
        possible_strings <- c(paste(w1, w2), paste(w2, w1))
        print(possible_strings)
        print(seen)
        curr_p <- all_p[seen == possible_strings[1]]
        if (length(curr_p) == 0) {
          curr_p <- all_p[seen == possible_strings[2]]
          
        }
        p_mat[w1, w2] <- curr_p
      }
    }
  }
  tree <- read.newick("phylogeny.nwk")
  pdf("new_version_of_Ines_manuscript/corrplot.pdf",
      height = 5, width = 10)
  tree_plot <- ggtree(tree, layout = "rectangular")  +
    scale_y_continuous(expand=expand_scale(0.17))  +
    theme(plot.margin = unit(c(-0.22,-2.3,2.2,0), "cm"))
  corr_plot <- ggcorrplot(cor_matrix, hc.order = FALSE)  +
    theme(plot.margin = unit(c(0,0,0,-2.3), "cm"),
            legend.text = element_text(size=10),
          legend.key.size = unit(1.5, 'cm'))
  for (w1 in 1:length(wasps_temp)) {
    for (w2 in 1:length(wasps_temp)) {
      corr_plot <- corr_plot + annotate("text", x = w1, y = w2, label = p_mat[w1, w2])
    }
  }
  grid.arrange(tree_plot, corr_plot, ncol = 2,
               widths = c(0.1,3))
  dev.off()
  par(mfrow = c(2,5))
  done <- c()
  for (species in 1:length(wasps)) {
    for (species2 in 1:length(wasps)) {
      if (!species == species2) {
        curr_string <- sort(c(wasps[species], wasps[species2]))
        curr_string <- paste(curr_string, collapse = "_")
        if (!curr_string %in% done) {
          done <- c(done, curr_string)
          corr <- cor.test(depths[, species2], depths[, species],
                           method = "spearman")
          plot(depths[, species2] ~ depths[, species],
               pch = "", xlab = wasps[species],
               ylab = wasps[species2],
               main = paste("rho = ", round(corr$estimate, 3),
                            "; P = ", signif(corr$p.value, 3),
                            sep = ""),
               cex.lab = 1.5)
          cols <- rep("black", length(depths[, species]))
          if (length(highlight > 0)) {
            cols[depth_ordered %in% highlight] <- "red"
          }
          text(x = depths[, species], y = depths[, species2],
               labels = circles_clean, col = cols)
        }
      }
    }
  }
}

# make a ridgeplot with the posteriors for all the regression coefficients
plot_posteriors <- function(data, focus, model, outcome, label, baseline, col_by_HIM = FALSE, subtract_for_HIM = FALSE, file_name = NA, pretty_names = NA, cube_root = FALSE) {
  # focus: the variable of interest in this comparison
  if (col_by_HIM == TRUE) {
    data$HIM[data$HIM == "y"] <- "HIM"
    data$HIM[data$HIM == "n"] <- "no HIM"
  }
  # extract the MCMC draws
  draws <- as.matrix(model)
  intercepts <- draws[,1]
  # if you want to combine the coefficients for the predictor of interest
  # with those for the HIM effect
  if (subtract_for_HIM == TRUE) {
    draws_HIM <- draws[, "b_HIMy"]
  }
  # get only the draws corresponding to
  # the coefficients for the predictor of interest
  # and clean up the parameter names
  draws <- draws[,str_starts(colnames(draws), paste("b_", focus, sep = ""))]
  categories <- substr(colnames(draws), (nchar(focus) + 3), nchar(colnames(draws)))
  if (length(pretty_names) > 1) {
    categories <- pretty_names
  }
  # make a vector with all the draws for the coefficients
  all_values <- c()
  all_categories <- c()
  # store the HIM status just for coloring in the distributions
  HIM <- c()
  for (category in categories) {
    values <- as.vector(draws[, categories == category])
    if (col_by_HIM == TRUE) {
      curr_him <- unique(data[data[,focus] == category, "HIM"])
      HIM <- c(HIM, rep(curr_him, length(values)))
      # if you want to combine with the HIM coefficients
      # and if the current segment has no HIM,
      # subtract the coefficient for HIM
      if (subtract_for_HIM == TRUE) {
        if (curr_him == "no HIM") {
          values <- values - draws_HIM
        }
      } 
    }
    all_values <- c(all_values, values)
    all_categories <- c(all_categories, rep(category, length(values)))
    all_categories <- gsub("susceptibilityn", "non-suitable host", all_categories)
  }
  # get differences for all possible pairwise comparisons between categories
  cat_names_all <- c(baseline, categories)
  for (cat in cat_names_all) {
    for (cat2 in cat_names_all) {
      # for the baseline, the prediction comes just from the intercept
      if (cat == baseline) {
        values1 <- exp(intercepts)
        if (cube_root == TRUE) {
          values1 <- intercepts ^ 3
        }
      }
      else {
        values1 <- exp(intercepts + draws[, categories == cat])
        if (cube_root == TRUE) {
          values1 <- (intercepts + draws[, categories == cat]) ^ 3
        }
      }
      if (cat2 == baseline) {
        values2 <- exp(intercepts)
        if (cube_root == TRUE) {
          values2 <- intercepts ^ 3
        }
      }
      else {
        values2 <- exp(intercepts + draws[, categories == cat2])
        if (cube_root == TRUE) {
          values2 <- (intercepts + draws[, categories == cat2]) ^ 3
        }
      }
      diffs <- values1 - values2
      # don't compare a category to itself
      if (!cat == cat2) {
        # probability that the mean is higher for category 1 than
        # category 2
        prob_above <- mean(diffs > 0)
        # plot out effect size only if the probability of an effect
        # in either direction is over 80%
        if (prob_above > 0.8 | prob_above < 0.2) {
          print("*********************")
          print(cat)
          print(cat2)
          print("Posterior median:")
          print(median(diffs))
          print(paste("95% HPDI for difference compared to", cat2))
          print(HPDinterval(as.mcmc(diffs), prob = 0.95))
          print("Prob > 0:")
          print(mean(diffs > 0)) 
        }
      }
    }
  }
  # change terminology from "segment" to "circle"
  if (focus == "Segment") {
    focus <- "Circle"
    all_categories <- gsub("Segment_", "Circle ", all_categories)
  }
  all_categories <- gsub("_", " ", all_categories)
  all_categories <- gsub("\\.", ". ", all_categories)
  data_for_plot <- data.frame(Outcome = all_values,
                              Predictor = all_categories)
  colnames(data_for_plot) <- c("Outcome", focus)
  # colour in the categories by HIM presence
  # only meaningful when the predictor of interest is
  # the circle
  if (col_by_HIM == TRUE) {
    data_for_plot$HIM <- HIM
    if (length(unique(HIM)) == 1) {
      colours <- c("pink3")
    }
    else {
      colours <- c("lightblue2", "pink3")
    }
  }
  # order categories by mean effect
  averages <- aggregate(data_for_plot[, "Outcome"] ~ data_for_plot[, focus], FUN = mean)
  ordered <- averages[order(averages[,2]),1]
  ordered[ordered == "C. sesamiaemombasa"] <- "C. sesamiae mombasa"
  data_for_plot[data_for_plot[,focus] == "C. sesamiaemombasa", focus] <- "C. sesamiae mombasa"
  data_for_plot[, focus] <- factor(data_for_plot[, focus], levels = ordered)
  focus <- sym(focus)
  # make ridge plot
  red <- rgb(186, 0, 0, maxColorValue = 255)
  blue <- rgb(115, 186, 215, maxColorValue = 255)
  if (col_by_HIM == TRUE) {
    ridgeplot <- ggplot(data_for_plot, aes(x = Outcome, y = !!focus, fill = HIM)) +
      geom_density_ridges_gradient(show.legend = TRUE) +
      theme_bw() + theme(legend.position = "right", legend.text = element_text(size=20)) +
      scale_fill_manual(guide = "legend", values=c("HIM"="pink3", "no HIM"="lightblue")) + 
      geom_vline(xintercept = 0) +
      guides(fill=guide_legend(title="")) + labs(x=outcome) +
      theme(axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            axis.title.y = element_blank())
  }
  else {
    ridgeplot <- ggplot(data_for_plot, aes(x = Outcome, y = !!focus, fill = after_stat(x))) +
      geom_density_ridges_gradient() +
      scale_fill_gradient(low = blue, high = red) + theme_bw() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0) + labs(x=outcome) +
      theme(axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            axis.title.y = element_blank())
  }
  if (length(categories) < 10) {
    pdf(file_name, width = 13, height = 7)
    ridgeplot <- ridgeplot +
      theme(axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            axis.title.y = element_blank())
    if (length(categories) == 5) {
      var_lab_mkd = c("non-suitable host")
      for (i in 1:4) {
        var_lab_mkd <- c(var_lab_mkd, paste("*", categories[i], "*", sep = ""))
      }
      ridgeplot <- ridgeplot + scale_y_discrete(labels = var_lab_mkd) + 
        theme(axis.text.y.left = element_markdown())
    } else {
      ridgeplot <- ridgeplot + theme(axis.text.y=element_text(face="italic"))
    }
  } else {
    pdf(file_name, width = 7, height = 10)
  }
  print(ridgeplot)
  dev.off()
  return(averages)
}

# make a ridgeplot with the posteriors for all the regression coefficients
plot_posteriors_difference <- function(data, focus, model, categories, pretty_names, title) {
  # set colours
  trans_red <- rgb(254, 132, 132, maxColorValue = 255, alpha = 10)
  trans_blue <- rgb(135, 206, 235, maxColorValue = 255, alpha = 10)
  red <- rgb(186, 0, 0, maxColorValue = 255)
  blue <- rgb(115, 186, 215, maxColorValue = 255)
  # get posterior means
  new_data <- data[1:length(categories),]
  new_data$hostwasp <- "B. fusca-C. sesamiae kitale"
  new_data$Segment <- "Segment_35"
  new_data$depth_cube_root_cpm_scaled <- 0
  new_data[, focus] <- categories
  set.seed(99)
  ped <- posterior_epred(model, newdata = new_data)
  ped_small <- ped[sample(1:nrow(ped), 1000),]
  line_colours <- rep(trans_blue, nrow(ped_small))
  line_colours[ped_small[,2] > ped_small[,1]] <- trans_red
  par("xpd" = TRUE)
  vioplot(log10(ped_small), col = c("white", "white"),
          xaxt = "n", yaxt = "n", ylim = c(-0.5,2.6),
          main = title, cex.main = 1.3)
  ticks <- c(0, 10, 100)
  axis(side = 2, at = log10(ticks + 1), labels = ticks, las = 2,
       cex.axis = 1.3)
  mtext("Est. mean # of integrations", side = 2, line = 3,
        cex = 1.3)
  text(x = c(0.8,2.2), y = -1.3, labels = pretty_names, srt = 15,
       cex = 1.3)
  for (draw in 1:nrow(ped_small)) {
    lines(log10(ped_small[draw,]) ~ c(1,2), col = line_colours[draw])
  }
  par("xpd" = FALSE)
}

# plot expected means per categorical predictor + depth
plot_posteriors_w_depth <- function(model, variable, baseline, file_name, data, height, selected, selected_pretty, maximum) {
  # generate some CPMs for the x-axis
  corr_depths <- seq(0, 1.3, by = 0.1)
  depth_cube_root <- corr_depths ^ (1/3)
  scale_mean <- mean(data$depth_cube_root_cpm)
  scale_sd <- sd(data$depth_cube_root_cpm)
  depth_cube_root_cpm_scaled <- (depth_cube_root - scale_mean)/scale_sd
  # generate predictions from model
  new_data <- data[1:(length(selected) * length(depth_cube_root)),]
  new_data$hostwasp <- "B. fusca-C. sesamiae kitale"
  new_data$Segment <- "Segment_35"
  new_data[, variable] <- rep(selected, each = length(depth_cube_root))
  new_data$depth_cube_root_cpm_scaled <- rep(depth_cube_root_cpm_scaled, 2)
  ped <- posterior_epred(model, newdata = new_data, re_formula = NA)
  pdf(file_name, width = 4, height = 4)
  plot(seq(0, maximum, length.out = length(depth_cube_root)) ~ depth_cube_root,
       pch = "", main = "", xlab = "Corrected DPM", 
       ylab = "Est. mean # of integrations",
       las = 2, xaxt = "n")
  labels <- c(0, 0.01, 0.05, 0.2, 0.5, 1)
  labels_loc <- labels ^ (1/3)
  axis(side = 1, at = labels_loc, labels = labels,
       las = 2)
  trans_red <- rgb(254, 132, 132, maxColorValue = 255, alpha = 10)
  trans_blue <- rgb(135, 206, 235, maxColorValue = 255, alpha = 10)
  red <- rgb(186, 0, 0, maxColorValue = 255)
  blue <- rgb(115, 186, 215, maxColorValue = 255)
  cols = c(trans_red, trans_blue)
  cols_dark = c(red, blue)
  selected_iter <- sample(1:nrow(ped), 1000)
  for (cat in 1:length(selected)) {
    for (iter in selected_iter) {
      lines(ped[iter, new_data[, variable] == selected[cat]] ~ depth_cube_root,
            col = cols[cat])
    }
    pred_means_across_draws <- apply(ped[selected_iter, new_data[, variable] == selected[cat]],
                                     MAR = 2, FUN = median)
    lines(pred_means_across_draws ~ depth_cube_root,
          col = cols_dark[cat], lwd = 2)
  }
  legend("topleft", legend = selected_pretty, col = cols_dark, lwd = 2,
         bty = "n")
  dev.off()  
  return()
}

# visualise posterior predictions as a function of a numerical variable of interest
plot_ppd <- function(ppd, draw_number, data, outcome, predictor, xlab, ylab, logged = FALSE, plot_medians = FALSE) {
  # sample the desired number of draws from the PPD
  ppd_plotting <- ppd[sample(1:dim(ppd)[1], draw_number, replace = FALSE),]
  # if the values should be log-transformed prior to plotting
  # (for the model of the number of integrations)
  if (logged == TRUE) {
    options(scipen = 1)
    # make a big vector with values for all the MCMC draws combined
    all_outcome <- c()
    all_pred <- c()
    for (draw in 1:dim(ppd_plotting)[1]) {
      all_outcome <- c(all_outcome, log10(ppd_plotting[draw,] + 1))
      all_pred <- c(all_pred, data[, predictor])
    }
    smoothScatter(all_outcome ~ all_pred, nrpoints = 0,
                  nbin = 500, ylab = ylab,
                  xlab = xlab,
                  yaxt = "n", xaxt = "n")
    ticks <- log10(c(1, 11, 101, 1001))
    axis(side = 2, at = ticks, labels = c(0, 10, 100, 1000),
         las = 2)
    labels <- c(0, 0.01, 0.05, 0.2, 0.5, 1)
    axis(side = 1, at = labels ^ (1/3), labels = labels,
         las = 2)
    points(log10(data[, outcome] + 1) ~ data[, predictor],
           pch = 20, col = "black")  
  }
  # if the values should not be log-transformed prior to plotting
  # (for the model of the number of integrations)
  else {
    all_outcome <- c()
    all_pred <- c()
    for (draw in 1:dim(ppd_plotting)[1]) {
      all_outcome <- c(all_outcome, ppd_plotting[draw,])
      all_pred <- c(all_pred, data[, predictor])
    }
    if (plot_medians == TRUE) {
      all_pred_round <- round(all_pred, digits = 1)
      simulated <- aggregate(all_outcome ~ all_pred_round, FUN = mean)
      real <- aggregate(data$depth_cube_root_cpm ~ round(data$log.host_BUSCO_depth, digits = 1), FUN = mean)
      plot(real[,2] ~ real[,1],
           type = "l", col = "steelblue3", lwd = 2, xlab = "log host BUSCO depth",
           ylab = "Corrected CPM", yaxt = "n")
      lines(simulated[,2] ~ simulated[,1],
            col = "black", lwd = 2)
      labels <- c(0, 10, 100, 300)
      axis(side = 2, at = labels ^ (1/3), labels = labels,
           las = 2)
      legend("bottomright", legend = c("Observed data", "Posterior predictions"),
             lwd = 2, col = c("black", "steelblue3"), bty = "n")
    } else {
      smoothScatter(log10(all_outcome) ~ all_pred, nrpoints = 0,
                    nbin = 500, ylab = ylab,
                    xlab = xlab,
                    las = 2, xlim = c(0, max(data[, predictor])),
                    yaxt = "n")
      labels <- c(0, 10, 100)
      axis(side = 2, at = log10(labels + 1), labels = labels,
      las = 2)
      points(log10(data[, outcome]) ~ data[, predictor],
             pch = 20, col = "black")  
    }
  }
}

# visualise posterior predictions as a function of segment ID
plot_ppd_by_segment <- function(ppd, draw_number, data, outcome, ylab, logged = FALSE) {
  # see plot_ppd for the structure of the code
  ppd_plotting <- ppd[sample(1:dim(ppd)[1], draw_number, replace = FALSE),]
  desert <- rgb(1, 0.89, 0.77, alpha = 0.5)
  trans_blue <- rgb(135, 206, 235, maxColorValue = 255, alpha = 100)
  if (logged == TRUE) {
    maximum_ppd <- max(log(ppd_plotting + 1)) + 0.2
    maximum_real <- max(log(data[, outcome] + 1)) + 0.2
    stripchart(log(ppd_plotting[1,] + 1) ~ data[, "Segment"],
               ylab = ylab,
               xlab = "", ylim = c(0, max(c(maximum_ppd, maximum_real))),
               pch = "", vertical = TRUE, 
               las = 2, cex.axis = 0.7, method = "jitter")
    values_to_plot <- c()
    segment_ids <- c()
    for (draw in 1:dim(ppd_plotting)[1]) {
      values_to_plot <- c(values_to_plot, log(ppd_plotting[draw,] + 1))
      segment_ids <- c(segment_ids, data[, "Segment"])
    }
    vioplot(values_to_plot ~ segment_ids,
            col = desert, 
            add = TRUE,
            rectCol = "grey",
            lineCol = "grey",
            border = "grey")
    stripchart(log(data[, outcome] + 1) ~ data[, "Segment"],
               pch = 20, col = trans_blue, vertical = TRUE, 
               add = TRUE, method = "jitter")  
  }
  else {
    par("mar" = c(5.5, 4.1, 4.1, 2.1))
    maximum_ppd <- max(ppd_plotting) + 0.2
    maximum_real <- max(data[, outcome]) + 0.2
    minimum_ppd <- min(ppd_plotting) - 0.2
    circle_names <- data[,"Segment"]
    circle_names <- gsub("Segment_", "Circle ", circle_names)
    stripchart(ppd_plotting[1,] ~ circle_names,
               ylab = ylab,
               xlab = "", ylim = c(minimum_ppd, max(c(maximum_ppd, maximum_real))),
               pch = "", vertical = TRUE, cex.axis = 0.7, las = 2)
    values_to_plot <- c()
    segment_ids <- c()
    for (draw in 1:dim(ppd_plotting)[1]) {
      values_to_plot <- c(values_to_plot, ppd_plotting[draw,])
      segment_ids <- c(segment_ids, data[, "Segment"])
    }
    vioplot(values_to_plot ~ segment_ids,
            col = desert, 
            add = TRUE,
            rectCol = "grey",
            lineCol = "grey",
            border = "grey")
    stripchart(data[, outcome] ~ data[, "Segment"],
               pch = 20, col = trans_blue, vertical = TRUE, add = TRUE) 
    medians <- aggregate(data[, outcome] ~ data[, "Segment"], FUN = median)
    points(medians, pch = 20, col = "steelblue")
  }
}

# plot the true value of a (simulated) parameter,
# as well as the 95% HPDI from brms
plot_true_with_estimated <- function(one_sim, positions, cols, labels, with_depth_only = FALSE,
                                     xlab = "Parameter value") {
  lows <- one_sim$HPDI_depth[1, positions]
  highs <- one_sim$HPDI_depth[2, positions]
  positions_no_depth <- positions
  positions_no_depth <- positions_no_depth - 1
  lows_no_depth <- one_sim$HPDI_no_depth[1, positions_no_depth]
  highs_no_depth <- one_sim$HPDI_no_depth[2, positions_no_depth]
  trues <- one_sim$real_values[positions]
  y <- seq(1, length(positions))
  minimum <- min(highs, lows, trues, highs_no_depth, lows_no_depth)
  maximum <- max(highs, lows, trues, highs_no_depth, lows_no_depth)
  if (with_depth_only == TRUE) {
    default <- par("mar")
    par("mar" = c(5.1, 7.1, 4.1, 2.1))
    par("xpd" = TRUE)
    plot(y ~ highs, pch = "", xlim = c(minimum, maximum),
         main = "", yaxt = "n",
         ylab = "",
         xlab = xlab, cex.lab = 0.75, 
         cex.axis = 0.75)
    text(x = -2.8, y = y,
         labels = labels, cex = 0.75)
    par("xpd" = FALSE)
    abline(v = 0, lty = 3, col = "grey", lwd = 2)
    par("mar" = default)
    segments(x0 = lows, y0 = y, x1 = highs, y1 = y)
    points(y ~ trues, col = cols, pch = 20)    
  } 
  else {    
    default <- par("mar")
    par("mar" = c(5.1, 7.1, 4.1, 2.1))
    par("xpd" = TRUE)
    plot(y ~ highs_no_depth, pch = "", xlim = c(minimum, maximum),
         main = "", yaxt = "n",
         ylab = "",
         xlab = xlab, cex.lab = 0.75, 
         cex.axis = 0.75)
    text(x = -2.8, y = y,
         labels = labels, cex = 0.75)
    par("xpd" = FALSE)
    abline(v = 0, lty = 3, col = "grey", lwd = 2)
    par("mar" = default)
    segments(x0 = lows_no_depth, y0 = y, x1 = highs_no_depth, y1 = y)
    points(y ~ trues, col = cols, pch = 20)    
  }
}

# given a categorical variable, set coefficients to 0 for all
# except for predefined categories, which will get a pre-defined coefficient
simulate_coefficients <- function(variable, to_modify, value_to_set_to, int_data_sim, baseline) {
  values <- levels(int_data_sim[, variable])
  coefs <- rep(0, dim(int_data_sim)[1])
  coefs_unique <- rep(0, length(values))
  for (cat in to_modify) {
    coefs[int_data_sim[, variable] == values[cat]] <- value_to_set_to
    coefs_unique[cat] <- value_to_set_to
  }
  coefs_unique <- coefs_unique[!values == baseline]
  return(list("coefs" = coefs, "coefs_unique" = coefs_unique))
} 

# make a simulated data set and see if depth is controlled for appropriately
simulate_data <- function(indata, seed, all_visuals = FALSE, weak_prior = FALSE) {
  set.seed(seed)
  # make sure the dodgy sample has been removed
  indata <- indata[indata$sample != "FAW6",]
  # based on posterior median in negative binomial model
  re_sigma <- 0.79
  int_data_sim <- indata
  sample_names <- levels(int_data_sim$sample)
  log_cpm <- log(indata$cpm)
  indata$log_cpm <- log_cpm
  # pick a value to use as the baseline depth
  mean_depth <- mean(indata$log_cpm[indata$cpm != 0])
  means <- rep(mean_depth, dim(int_data_sim)[1])
  # increase or decrease depth for certain Segments/hosts/wasps
  coefs <- simulate_coefficients("Segment", c(4, 7, 10), 0.5, int_data_sim, "Segment_35")[["coefs"]]
  means <- means + coefs
  coefs <- simulate_coefficients("hostwasp", c(4), 0.5, int_data_sim, "B. fusca-C. sesamiae kitale")[["coefs"]]
  means <- means + coefs
  # sample random values from around these means
  std <- sd(indata$log_cpm[indata$cpm != 0])
  sim_depths <- rnorm(dim(int_data_sim)[1], mean = means, sd = std)
  int_data_sim$cpm <- exp(sim_depths)
  # take cube root
  int_data_sim$depth_cube_root_cpm <- int_data_sim$cpm ^ (1/3)
  # check that simulated depths look reasonable and show
  # the specified effects
  if (all_visuals == TRUE) {
    default <- par("mar")
    par(mfrow = c(1,1))
    par("mar" = c(7.1, 4.1, 4.1, 2.1))
    vioplot(depth_cube_root_cpm ~ hostwasp, cex.axis = 0.7, col = "lightyellow",
            ylab = "Depth", xlab = "", las = 2, data = int_data_sim)
    vioplot(depth_cube_root_cpm ~ Segment, cex.axis = 0.7, col = "lightyellow",
            ylab = "Depth", xlab = "", las = 2, data = int_data_sim)
    par("mar" = c(9.1, 4.1, 4.1, 2.1))
    par("mar" = default)  
    plot(density(log(indata[indata$host == "C. partellus",]$depth_cube_root_cpm)), main = "", xlab = "ln Depth",
         col = "skyblue",
         lwd = 2)
    lines(density(log(indata$depth_cube_root_cpm)), col = "steelblue", lwd = 2)
    legend("topright", legend = c("observed", "simulated"), col = c("skyblue", "steelblue"),
           lwd = 2)
  }
  # model the number of integrations as a function of depth
  # simulate new numbers of integrations
  model_int <- glm.nb(int ~ depth_cube_root_cpm, data = indata)
  intercept <- coef(model_int)[1]
  means <- intercept + (int_data_sim$depth_cube_root_cpm * coef(model_int)[2])
  betas <- c(coef(model_int)[2])
  # shift the number of integrations based on predictor values
  # (i.e. simulate effects)
  coefs_segment <- simulate_coefficients("Segment", c(2,3), 1, int_data_sim, "Segment_35")
  betas <- c(betas, coefs_segment[["coefs_unique"]])
  means <- means + coefs_segment[["coefs"]]
  coefs_host <- simulate_coefficients("hostwasp", c(5), 1, int_data_sim, "B. fusca-C. sesamiae kitale")
  betas <- c(betas, coefs_host[["coefs_unique"]])
  means <- means + coefs_host[["coefs"]]
  # simulate sample-level random effects
  alphas <- rnorm(length(sample_names), mean = 0, sd = re_sigma)
  counter <- 1
  for (s in sample_names) {
    means[int_data_sim$sample == s] <- means[int_data_sim$sample == s] + alphas[counter]
    counter <- counter + 1
  }
  # simulate numbers of integrations
  phi <- model_int$theta
  int_data_sim$int <- rnbinom(dim(int_data_sim)[1], 
                              mu = exp(means),
                              size = phi)
  # run model either controlling for depth or not controlling for it
  prior_b0 <- set_prior("normal(0,3)", class = "Intercept")
  prior_b <- set_prior("normal(0,2)", class = "b")
  prior_sd <- set_prior("normal(0,5)", class = "sd")
  prior_shape <- set_prior("inv_gamma(0.4, 0.3)", class = "shape")
  priors <- c(prior_b0, prior_b, prior_sd, prior_shape)
  depth_model <- brm(int ~ depth_cube_root_cpm + Segment + hostwasp + (1|sample),
                   data = int_data_sim, family = negbinomial(),
                   prior = priors, cores = 4, iter = 500, 
                   warmup = 250, silent = 2, refresh = 0)
  no_depth_model <- brm(int ~ Segment + hostwasp + (1|sample),
                     data = int_data_sim, family = negbinomial(),
                     prior = priors, cores = 4, iter = 500, 
                     warmup = 250, silent = 2, refresh = 0)
  if (all_visuals == TRUE) {
    print(summary(depth_model))
    par(mfrow = c(1,2))
    plot(log(int + 1) ~ depth_cube_root_cpm, data = indata,
         main = "Observed")
    plot(log(int + 1) ~ depth_cube_root_cpm, data = int_data_sim,
         main = "Simulated") 
    print(cor.test(log(indata$int + 1), indata$depth_cube_root_cpm))
    print(cor.test(log(int_data_sim$int + 1), int_data_sim$depth_cube_root_cpm))
  }
  # store various statistics
  draws_no_depth <- as.matrix(no_depth_model)
  draws_depth <- as.matrix(depth_model)
  stat_names <- c(colnames(draws_depth)[startsWith(colnames(draws_depth), "b_")],
                  "shape",
                  colnames(draws_depth)[startsWith(colnames(draws_depth), "r_sample")],
                  "sd_sample__Intercept")
  stat_names_no_depth <- c(colnames(draws_no_depth)[startsWith(colnames(draws_no_depth), "b_")],
                           "shape",
                           colnames(draws_no_depth)[startsWith(colnames(draws_no_depth), "r_sample")],
                           "sd_sample__Intercept")
  real_values <- c(intercept,
                   betas,
                   phi,
                   alphas,
                   re_sigma)
  posterior_medians_no_depth <- apply(draws_no_depth[, stat_names_no_depth], 2, median)
  posterior_medians_depth <- apply(draws_depth[, stat_names], 2, median)
  HPDI_no_depth <- apply(draws_no_depth[, stat_names_no_depth], 2, FUN = function(x) {HPDinterval(as.mcmc(x), prob = 0.95)})
  HPDI_depth <- apply(draws_depth[, stat_names], 2, FUN = function(x) {HPDinterval(as.mcmc(x), prob = 0.95)})
  above0 <- apply(draws_depth[, stat_names], 2, FUN = function(x) {mean(x > 0)})
  above0_no_depth <- apply(draws_no_depth[, stat_names_no_depth], 2, FUN = function(x) {mean(x > 0)})
  if (all_visuals == TRUE) {
    for (param in 1:length(stat_names)) {
      print("*****")
      print(stat_names[param])
      print(real_values[param])
      print(HPDI_depth[, param])
      print(above0[param])
    }
  }
  return(list("stat_names" = stat_names,
              "stat_names_no_depth" = stat_names_no_depth,
              "real_values" = real_values,
              "posterior_medians_no_depth" = posterior_medians_no_depth,
              "posterior_medians_depth" = posterior_medians_depth,
              "HPDI_no_depth" = HPDI_no_depth,
              "HPDI_depth" = HPDI_depth,
              "above0" = above0,
              "above0_no_depth" = above0_no_depth))
}

# wrapper around simulate data
# performs sim_number simulations
# for each simulation, checks for each parameter whether the true value
# is within the 95% HPDI
simulate_data_wrapper <- function(int_data, sim_number, seed, weak_prior = FALSE) {
  # create an empty matrix, with simulations in rows and parameters in
  # columns
  # make two versions of it, one for the model that controls for depth
  # and one for the model that doesn't
  captured <- matrix(0L, nrow = sim_number, ncol = 46)
  captured_no_depth <- matrix(0L, nrow = sim_number, ncol = 46)
  # create matrices to store the probability of each estimate 
  # being above 0 for each parameter
  # this is useful to check that depth is controlled for successfully
  above0 <- matrix(0L, nrow = sim_number, ncol = 46)
  above0_no_depth <- matrix(0L, nrow = sim_number, ncol = 46)
  # and posterior medians
  post_medians <- matrix(0L, nrow = sim_number, ncol = 46)
  post_medians_no_depth <- matrix(0L, nrow = sim_number, ncol = 46)
  for (sim in 1:sim_number) {
    print("Simulation:")
    print(sim)
    # perform simulation
    options(warn = -1)
    out <- simulate_data(int_data, seed + sim, FALSE, weak_prior = weak_prior)
    options(warn = 0)
    # check if the sad wasps were detected
    # fill in the corresponding rows of the output matrices
    for (param in 1:46) {
      captured[sim, param] <- out$real_values[param] > out$HPDI_depth[1, param] & out$real_values[param] < out$HPDI_depth[2, param] 
      above0[sim, param] <- out$above0[param]
      post_medians[sim, param] <- out$posterior_medians_depth[param]
      if (!param == 2) {
        if (param == 1) {
          captured_no_depth[sim, param] <- out$real_values[param] > out$HPDI_no_depth[1, param] & out$real_values[param] < out$HPDI_no_depth[2, param] 
          above0_no_depth[sim, param] <- out$above0_no_depth[param]
          post_medians_no_depth[sim, param] <- out$posterior_medians_no_depth[param]
        }
        else {
          captured_no_depth[sim, param] <- out$real_values[param] > out$HPDI_no_depth[1, (param - 1)] & out$real_values[param] < out$HPDI_no_depth[2, (param - 1)] 
          above0_no_depth[sim, param] <- out$above0_no_depth[param - 1]
          post_medians_no_depth[sim, param] <- out$posterior_medians_no_depth[param - 1]
        }
      }
    }
  }
  return(list("captured" = captured, "captured_no_depth" = captured_no_depth,
              "above0" = above0, "above0_no_depth" = above0_no_depth,
              "post_medians" = post_medians, 
              "post_medians_no_depth" = post_medians_no_depth))
}

######

# Read in and clean data
######
int_data <- read.csv("chimera_depth_Bayes_18_03_25.txt", sep = "\t")
ls_data <- read.csv("depth_nb_reads_30.05.txt", sep = "\t")
# extract library size data and tack it on as an extra column in
# the main data frame
ls <- c()
for (line in 1:dim(int_data)[1]) {
  ls <- c(ls, ls_data$nb_reads[ls_data$sample == int_data[line, 2]])
}
int_data$library_size <- ls
# initial file had incorrect library size
int_data$library_size[int_data$sample == "a1"] <- 352004484

# clean data

# convert the number of integrations to an integer
# for modelling as a count
int_data$int <- as.integer(round(int_data$int))
# convert categorical predictors to factors
int_data[int_data$Segment == "Segment_20/33", "Segment"] <- "Segment_20_33"
segment_names <- unique(int_data$Segment)
int_data$Segment = factor(int_data$Segment, levels = segment_names)
# give prettier names to the wasps
original_wasps <- unique(int_data$wasp)
new <- c("C. sesamiae kitale",
         "C. sesamiae mombasa",
         "C. flavipes",
         "C. icipe",
         "C. typhae")
for (wasp in 1:length(original_wasps)) {
  int_data$wasp[int_data$wasp == original_wasps[wasp]] <- new[wasp]
}
int_data$wasp = factor(int_data$wasp, levels = new)
# give prettier names to hosts
original <- unique(int_data$host)
new <- c("B. fusca",
         "S. calamistis",
         "C. partellus",
         "S. frugiperda",
         "S. nonagrioides")
for (host in 1:length(original)) {
  int_data$host[int_data$host == original[host]] <- new[host]
}
int_data$host = factor(int_data$host, levels = new)

int_data$sample <- as.factor(int_data$sample)
# set aside a version of the data frame that has both circles with HIMs and without
int_data_w_HIM <- int_data
# only keep circles that have an HIM in the main data frame
int_data <- int_data[int_data$HIM == "y",]
int_data$Segment <- droplevels(int_data$Segment)

# subtract wasp depth from the total depth to approximate
# the depth for host alone
int_data$norm_depth1 <- int_data$depth - int_data$wasp_BUSCO_depth
# set negative values to 0
int_data$norm_depth1[int_data$norm_depth1 < 0] <- 0
# do the same for the data frame that has circles with no HIM
int_data_w_HIM$norm_depth1 <- int_data_w_HIM$depth - int_data_w_HIM$wasp_BUSCO_depth
int_data_w_HIM$norm_depth1[int_data_w_HIM$norm_depth1 < 0] <- 0
# make a new variable that combines the host and wasp
# information
int_data$hostwasp <- paste(int_data$host, int_data$wasp, sep = "-")
int_data_w_HIM$hostwasp <- paste(int_data_w_HIM$host, int_data_w_HIM$wasp, sep = "-")

# making sure that the first level is always the correct
# baseline for factors
other_levels <- unique(as.character(int_data$Segment))
other_levels <- other_levels[other_levels != "Segment_35"]
int_data$Segment <- factor(int_data$Segment, 
                                 levels = c("Segment_35", other_levels))
other_levels <- unique(as.character(int_data$host))
other_levels <- other_levels[other_levels != "S. calamistis"]
int_data$host <- factor(int_data$host, levels = c("S. calamistis", 
                                                              other_levels))
other_levels <- unique(as.character(int_data$wasp))
other_levels <- other_levels[other_levels != "C. sesamiae kitale"]
int_data$wasp <- factor(int_data$wasp, levels = c("C. sesamiae kitale",
                                                              other_levels))
int_data$host_susceptibility <- factor(int_data$host_susceptibility,
                                             levels = c("y", "n"))
other_levels <- unique(as.character(int_data$hostwasp))
other_levels <- other_levels[other_levels != "B. fusca-C. sesamiae kitale"]
int_data$hostwasp <- factor(int_data$hostwasp, levels = c("B. fusca-C. sesamiae kitale",
                                                                          other_levels))
other_levels <- unique(as.character(int_data_w_HIM$Segment))
other_levels <- other_levels[other_levels != "Segment_35"]
int_data_w_HIM$Segment <- factor(int_data_w_HIM$Segment, 
                                 levels = c("Segment_35", other_levels))
other_levels <- unique(as.character(int_data_w_HIM$host))
other_levels <- other_levels[other_levels != "S. calamistis"]
int_data_w_HIM$host <- factor(int_data_w_HIM$host, levels = c("S. calamistis", 
                                                              other_levels))
other_levels <- unique(as.character(int_data_w_HIM$wasp))
other_levels <- other_levels[other_levels != "C. sesamiae kitale"]
int_data_w_HIM$wasp <- factor(int_data_w_HIM$wasp, levels = c("C. sesamiae kitale",
                                                              other_levels))
int_data_w_HIM$host_susceptibility <- factor(int_data_w_HIM$host_susceptibility,
                                             levels = c("y", "n"))
other_levels <- unique(as.character(int_data_w_HIM$hostwasp))
other_levels <- other_levels[other_levels != "B. fusca-C. sesamiae kitale"]
int_data_w_HIM$hostwasp <- factor(int_data_w_HIM$hostwasp, levels = c("B. fusca-C. sesamiae kitale",
                                                              other_levels))



# normalise for library size
int_data$cpm <- int_data$norm_depth1/(int_data$library_size/1000000)
int_data_w_HIM$cpm <- int_data_w_HIM$norm_depth1/(int_data_w_HIM$library_size/1000000)

# transform the normalised depth by taking the cube root
int_data$depth_cube_root <- int_data$norm_depth1 ^ (1/3)
int_data_w_HIM$depth_cube_root <- int_data_w_HIM$norm_depth1 ^ (1/3)
int_data$depth_cube_root_cpm <- int_data$cpm ^ (1/3)
int_data_w_HIM$depth_cube_root_cpm <- int_data_w_HIM$cpm ^ (1/3)
plot(depth_cube_root_cpm ~ depth_cube_root, data = int_data)

# visualise the effects of the library size normalisation
par(mfrow = c(1,2))
plot(log(int + 1) ~ depth_cube_root, data = int_data,
     ylab = "log(# int + 1)", xlab = "Corrected Depth ^ 1/3",
     main = "No LS normalisation")
plot(log(int + 1) ~ depth_cube_root_cpm, data = int_data,
     ylab = "log(# int + 1)", xlab = "CPM ^ 1/3",
     main = "LS normalisation")

# compare library sizes between samples
par(mfrow = c(1,1))
ls_table <- aggregate(log(library_size) ~ sample,
                      FUN = unique,
                      data = int_data)
barplot(ls_table[,2] ~ ls_table[,1])

# check which samples have which segments
for (s in unique(int_data_w_HIM$sample)) {
  print("***")
  print(s)
  print(length(unique(int_data_w_HIM[int_data_w_HIM$sample == s, "Segment"])))
  print(unique(int_data_w_HIM[int_data_w_HIM$sample == s, "Segment"]))
}

######

# Check whether the order of segments (in terms of depth/number of integrations)
# is conserved.
######
# for this analysis, only keep segments that are present in all the species
seg_number <- aggregate(depth_cube_root ~ Segment, data = int_data_w_HIM,
          FUN = length)
seg_to_keep <- seg_number[seg_number[,2] == 21, 1]
int_data_w_HIM_all_species <- int_data_w_HIM[int_data_w_HIM$Segment %in% seg_to_keep,]
# calculate mean depth per segment
depth_summaries <- aggregate(depth_cube_root_cpm ~ Segment, data = int_data_w_HIM_all_species,
                         FUN = mean)
# order the segments by depth
depth_ordered_temp <- depth_summaries[order(depth_summaries[,2], decreasing = FALSE),]
depth_ordered_temp <- cbind(depth_ordered_temp, seq(1:25))
colnames(depth_ordered_temp)[3] <- "Order"
depth_ordered <- depth_summaries[order(depth_summaries[,2], decreasing = FALSE), 1]
# make a vector of wasp IDs
wasps <- unique(int_data_w_HIM$wasp)

# plot correlations between depth values
plot_order_corrs(depth_ordered, int_data_w_HIM_all_species, wasps,
                 highlight = HIM_segments)
######

# Model the number of integrations 
######
# remove sample where oviposition appears to have 
# been unsuccessful
int_data <- int_data[!int_data$sample == "FAW6",]
int_data_w_HIM <- int_data_w_HIM[!int_data_w_HIM$sample == "FAW6",]
int_data <- droplevels(int_data)
int_data_w_HIM <- droplevels(int_data_w_HIM)

# run negative binomial model
prior_b0 <- set_prior("normal(0,3)", class = "Intercept")
prior_b <- set_prior("normal(0,2)", class = "b")
prior_sd <- set_prior("normal(0,5)", class = "sd")
prior_shape <- set_prior("inv_gamma(0.4, 0.3)", class = "shape")
priors <- c(prior_b0, prior_b, prior_sd, prior_shape)
int_data$depth_cube_root_cpm_scaled <- scale(int_data$depth_cube_root_cpm)
nb_model <- brm(int ~ Segment + hostwasp + depth_cube_root_cpm_scaled + (1|sample),
                 data = int_data, family = negbinomial(),
                 prior = priors, cores = 4, iter = 10000,
                warmup = 1000, seed = 90)

# difference between the number of integrations at norm depth Q1 vs Q3
summary(int_data$cpm)
depth_effect_HPDI(nb_model, depths = c(0.08108, 0.31570), 
                  int_data$depth_cube_root_cpm)

# the effect of depth shown with the effect of segment
ordered <- plot_posteriors_w_depth(nb_model, "Segment", "Segment_35", "new_version_of_Ines_manuscript/Figure_3b.pdf",
                                int_data, 5, selected = c("Segment_1", "Segment_14"),
                                selected_pretty = c("Circle 1", "Circle 14"), maximum = 150)
# effect of host suitability, shown with the effect of segment
plot_posteriors_w_depth(nb_model, "host_susceptibility", "y", "new_version_of_Ines_manuscript/Figure_5A.pdf",
                                        int_data, 5, selected = c("y", "n"),
                             selected_pretty = c("Suitable host", "Non-suitable host"), 
                             maximum = 100)

# check whether having the two non-suitable at the bottom is significant
empirical_dist <- permutations(n = length(levels(int_data$hostwasp)), r = length(levels(int_data$hostwasp)), v = levels(int_data$hostwasp))
nonsuitable <- as.character(unique(int_data$hostwasp[int_data$host_susceptibility == "n"]))
positive <- 0
for (row in 1:nrow(empirical_dist)) {
  if (empirical_dist[row, 7] %in% nonsuitable & empirical_dist[row, 8] %in% nonsuitable) {
    positive <- positive + 1
  }
}
positive/nrow(empirical_dist)

# obtain posterior predictive distribution and visualise either by depth 
# or by segment ID
set.seed(7)
ppd <- posterior_predict(nb_model)
pdf("new_version_of_Ines_manuscript/Figure_S7.pdf", width = 4, height = 4)
plot_ppd(ppd, 1000, int_data, "int", "depth_cube_root_cpm",
         "Corrected DPM", 
         "# of integrations", logged = FALSE)
dev.off()

# show results for circle, host and wasp
ordered <- plot_posteriors(int_data, "Segment", nb_model,
                           "Difference to Circle 35 in log # of integrations", "", "Segment_35", 
                           file_name = "new_version_of_Ines_manuscript/Figure_4.pdf")
x <- plot_posteriors(int_data, "hostwasp", nb_model,
                     expression("Difference to "~italic("C. s.")~" Kitale -"~italic("B. fusca")~" in log # of integrations"), "", "B. fusca-C. s. kitale", 
                     file_name = "new_version_of_Ines_manuscript/Figure_S9.pdf",
                     pretty_names = c("C. s. Mombasa - B.fusca",
                                      "C. s. Mombasa - S. calamistis",
                                      "C. s. Kitale - S. calamistis",
                                      "C. flavipes - C. partellus",
                                      "C. flavipes - S. frugiperda",
                                      "C. icipe - S. frugiperda",
                                      "C. typhae - S. nonagrioides"))

# modelling results
pdf("new_version_of_Ines_manuscript/Figure_5A.pdf", width = 6, height = 4)
par(mfrow = c(1,2))
plot_posteriors_difference(int_data, "hostwasp", nb_model, 
                           c("B. fusca-C. sesamiae mombasa", "S. calamistis-C. sesamiae mombasa"), 
                           pretty_names = c(expression(paste(italic("B. fusca"))), expression(paste(italic("S. calamistis")))),
                           title = substitute(paste(italic("C. s."), " Mombasa with...")))
plot_posteriors_difference(int_data, "hostwasp", nb_model, 
                           c("S. frugiperda-C. flavipes", "C. partellus-C. flavipes"), 
                           pretty_names = c(expression(paste(italic("S. frugiperda"))), expression(paste(italic("C. partellus")))),
                           title = substitute(paste(italic("C. flavipes"), " with...")))
dev.off()

# descriptive data analysis
# mombasa and flavipes when paired with a suitable vs
# non-suitable host
# set up colours
red <- rgb(186, 0, 0, maxColorValue = 255)
blue <- rgb(115, 186, 215, maxColorValue = 255)
pdf(file = "new_version_of_Ines_manuscript/Figure_5B.pdf", width = 6, height = 4)
par(mfrow = c(1,2))
# get data only for momabasa
mombasa <- int_data[int_data$wasp == "C. sesamiae mombasa",]
# get number of integrations per segment
# separately for suitable and non-suitable host
per_segment_data <- aggregate(int ~ Segment + host_susceptibility, data = mombasa,
                              FUN = median)
# make graph
yvalues <- c(0, 10, 60)
yvalues_logged <- log(yvalues + 1)
par("xpd" = TRUE)
suitability <- factor(per_segment_data$host_susceptibility, 
                         levels = c("n", "y"))
stripchart(per_segment_data[,3] ~ suitability, data = mombasa, vertical = TRUE,
           method = "jitter", main = substitute(paste(italic("C. s."), " Mombasa with...")),
           xlab = "",
           ylab = "# of integrations",
           col = c(blue, red),
           ylim = c(0,90),
           pch = 1, xaxt = "n", cex.main = 1.3, 
           cex.axis = 1.3, cex.lab = 1.3, las = 2)
text(x = c(0.8,1.8), y = -19, labels = c(expression(paste(italic("B. fusca"))), expression(paste(italic("S. calamistis")))),
     srt = 20, cex = 1.3)
# now go through all the same steps with flavipes
flavipes <- int_data[int_data$wasp == "C. flavipes",]
per_segment_data <- aggregate(int ~ Segment + host_susceptibility, data = flavipes,
                              FUN = median)
suitability <- factor(per_segment_data$host_susceptibility, 
                         levels = c("n", "y"))
stripchart(per_segment_data[,3] ~ suitability, data = mombasa, vertical = TRUE,
           method = "jitter", main = substitute(paste(italic("C. flavipes"), " with...")),
           xlab = "", pch = 1,
           ylab = "# of integrations",
           col = c(blue, red),
           ylim = c(0,90), xaxt = "n", cex.main = 1.3, 
           cex.axis = 1.3, cex.lab = 1.3, las = 2)
text(x = c(0.8,1.8), y = -19, labels = c(expression(paste(italic("S. frugiperda"))), expression(paste(italic("C. partellus")))),
     srt = 20, cex = 1.3)
par("xpd" = FALSE)
dev.off()

######

# Create simulated data sets and see if the model can recover the ground truth
######
# remove sample where oviposition appears to have 
# been unsuccessful
int_data <- int_data[!int_data$sample == "FAW6",]
int_data_w_HIM <- int_data_w_HIM[!int_data_w_HIM$sample == "FAW6",]
int_data <- droplevels(int_data)
int_data_w_HIM <- droplevels(int_data_w_HIM)

# do a single simulation to have access to the column names
one_sim <- simulate_data(int_data, 40, all_visuals = TRUE)

# do 100 simulations
captured_out <- simulate_data_wrapper(int_data, 100, 25)
# label the results from the 100 simulations
colnames(captured_out$captured) <- one_sim$stat_names
colnames(captured_out$captured_no_depth) <- one_sim$stat_names
colnames(captured_out$above0) <- one_sim$stat_names
colnames(captured_out$above0_no_depth) <- one_sim$stat_names
# write the results of the simulations to file
write.table(captured_out$captured, "captured_depth_120325.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(captured_out$captured_no_depth, "captured_no_depth_120325.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(captured_out$above0, "above0_depth_120325.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(captured_out$above0_no_depth, "above0_no_depth_120325.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
# read simulation results from file
captured <- read.table("captured_depth_120325.txt", header = TRUE)
captured_no_depth <- read.table("captured_no_depth_120325.txt", header = TRUE)
above0 <- read.table("above0_depth_120325.txt", header = TRUE)
above0_no_depth <- read.table("above0_no_depth_120325.txt", header = TRUE)

# for convenience, repeat the part of the simulation that creates the ground truth
# just to have the values easily available
depth_coefs_seg <- simulate_coefficients("Segment", c(4, 7, 10), 0.5, int_data, "Segment_35")[["coefs_unique"]]
depth_coefs_hostwasp <- simulate_coefficients("hostwasp", 4, 0.5, int_data, "B. fusca-C. sesamiae kitale")[["coefs_unique"]]
depth_coefs <- c(depth_coefs_seg, depth_coefs_hostwasp)
names(depth_coefs) <- one_sim$stat_names[startsWith(one_sim$stat_names, "b_")][-c(1,2)]

coefs_segment <- simulate_coefficients("Segment", c(2,3), 1, int_data, "Segment_35")[["coefs_unique"]]
coefs_hostwasp <- simulate_coefficients("hostwasp", c(5), 1, int_data, "B. fusca-C. sesamiae kitale")[["coefs_unique"]]
int_coefs <- c(coefs_segment, coefs_hostwasp)
names(int_coefs) <- one_sim$stat_names[startsWith(one_sim$stat_names, "b_")][-c(1,2)]
  
change_to_depth <- which(depth_coefs != 0)
change_to_int <- which(int_coefs != 0)

# plot a single simulation
beta_pos <- which(startsWith(colnames(above0), "b_"))[-c(1,2)]
pdf("new_version_of_Ines_manuscript/Figure_S11B.pdf", width = 4, height = 4)
par(mfrow = c(1,1))
pretty_names <- c("Circle 1", "Circle 10", expression(italic("C. flavipes-C. partellus")),
                  "Circle 11",
                  "Circle 16", "Circle 24",
                  expression(italic("C. s. kitale-S. calamistis")))
pretty_names_no_expr <- c("Circle 1", "Circle 10", "C. flavipes-C. partellus",
                          "Circle 11",
                          "Circle 16", "Circle 24",
                          "C. s. kitale-S. calamistis")
cols <- c(rep("steelblue3", 3), rep("grey", 4))
plot_true_with_estimated(one_sim, beta_pos[c(change_to_int, change_to_depth)], 
                         cols,
                         labels = pretty_names,
                         with_depth_only = TRUE,
                         xlab = "Coefficient")
dev.off()

# without controlling for depth
pdf("new_version_of_Ines_manuscript/Figure_S12B.pdf", width = 4, height = 4)
par(mfrow = c(1,1))
plot_true_with_estimated(one_sim, beta_pos[c(change_to_int, change_to_depth)], 
                         cols,
                         labels = pretty_names,
                         with_depth_only = FALSE,
                         xlab = "Coefficient")
dev.off()

# plot all simulations
pdf("new_version_of_Ines_manuscript/Figure_S11A.pdf", height = 4, width = 5)
columns_to_check <- beta_pos[c(change_to_int, change_to_depth)]
above <- apply(above0[, columns_to_check], 2, FUN = median) * 100
aboveq1 <- apply(above0[, columns_to_check], 2, FUN = function(x) {quantile(x, 1/3)}) * 100
aboveq3 <- apply(above0[, columns_to_check], 2, FUN = function(x) {quantile(x, 2/3)}) * 100
below <- -(100 - above)
belowq1 <- -(100 - aboveq1)
belowq3 <- -(100 - aboveq3)
for_plotting <- data.frame("Probability of effect" = c(above, below),
                           "Q1" = c(aboveq1, belowq1),
                           "Q3" = c(aboveq3, belowq3),
                           "Direction" = c(rep("Increased", length(above)),
                                           rep("Decreased", length(below))),
                           "Predictor" = rep(pretty_names_no_expr, 2),
                           check.names = FALSE)
for_plotting$Predictor <- factor(for_plotting$Predictor, levels = pretty_names_no_expr)
ggplot(for_plotting, aes(x = Predictor, group = Direction, y = `Probability of effect`, fill = Direction)) + 
  geom_bar(stat="identity", position= "identity") + coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c(Decreased = "skyblue", Increased = "salmon")) +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.title=element_blank()) +
  geom_hline(yintercept = c(95, -95), 
             color = "grey", linetype = "dashed", linewidth = 0.5) + 
  scale_x_discrete(labels = pretty_names)
dev.off()

pdf("new_version_of_Ines_manuscript/Figure_S12A.pdf", height = 4, width = 5)
above <- apply(above0_no_depth[, columns_to_check], 2, FUN = median) * 100
below <- -(100 - above)
for_plotting <- data.frame("Probability of effect" = c(above, below),
                           "Direction" = c(rep("Increased", length(above)),
                                           rep("Decreased", length(below))),
                           "Predictor" = rep(pretty_names_no_expr, 2),
                           check.names = FALSE)
for_plotting$Predictor <- factor(for_plotting$Predictor, levels = pretty_names_no_expr)
ggplot(for_plotting, aes(x = Predictor, group = Direction, y = `Probability of effect`, fill = Direction)) + 
  geom_bar(stat="identity", position= "identity") + coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c(Decreased = "skyblue", Increased = "salmon")) +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.title=element_blank()) +
  geom_hline(yintercept = c(95, -95), 
             color = "grey", linetype = "dashed", linewidth = 0.5) + 
  scale_x_discrete(labels = pretty_names)
dev.off()

pdf("new_version_of_Ines_manuscript/Figure_S10.pdf", height = 5, width = 8)
margins <- par("mar")
par("mar" = c(10.1, 4.1, 4.1, 2.1))
clean_names <- c("Intercept" , "Corrected CPM ^ (1/3)" , "Circle 1", "Circle 10" , "Circle 11" , "Circle 12" , "Circle 14" , "Circle 16" , "Circle 17" , "Circle 18" , "Circle 24" , "Circle 26" , "Circle 27" , "Circle 28" , "Circle 32" , "Circle 4", "Circle 7", "C. s. mombasa - B. fusca" , "C. s. mombasa - S. calamistis", "C. s. kitale - S. calamistis" , "C. flavipes - C. partellus", "C. flavipes - S. frugiperda" , "C. icipe - S. frugiperda", "C.typhae - S. nonagrioides" , "phi" , "alpha 1", "alpha 2" , "alpha 3" , "alpha 4" , "alpha 5" , "alpha 6" , "alpha 7" , "alpha 8" , "alpha 9" , "alpha 10" , "alpha 11" , "alpha 12" , "alpha 13" , "alpha 14", "alpha 15", "alpha 16" , "alpha 17" , "alpha 18" , "alpha 19" , "alpha 20" , "sigma")
barplot(colMeans(captured), 
        main = "",
        col = "salmon", names = clean_names,
        las = 2, cex.names = 0.75,
        cex.axis = 0.75, cex.lab = 0.75,
        ylab = "Proportion of simulations where 95% HPDI includes true value")
abline(h = 0.95, lty = 2, lwd = 2, col = "black")
par("mar" = margins)
dev.off()

# check how reliably the efects are detected
colnames(above0) <- one_sim$stat_names
median(above0$b_hostwaspC.partellusMC.flavipes)
median(above0$b_SegmentSegment_1)
median(above0$b_SegmentSegment_10)


######

# Model depth as outcome
######
# remove sample where oviposition appears to have 
# been unsuccessful
int_data <- int_data[!int_data$sample == "FAW6",]
int_data_w_HIM <- int_data_w_HIM[!int_data_w_HIM$sample == "FAW6",]
int_data <- droplevels(int_data)
int_data_w_HIM <- droplevels(int_data_w_HIM)

# set priors
prior_b0_depth <- set_prior("normal(0,4)",
                      class = "Intercept")
prior_b_depth <- set_prior("normal(0,2)",
                     class = "b")
prior_sd_depth <- set_prior("normal(0,2)", class = "sigma")
prior_re_sd_depth <- set_prior("normal(0,5)", class = "sd")

# fit model without the replication unit and HIM presence as a predictors
depth_model_no_HIM_pred <- brm(depth_cube_root_cpm ~ hostwasp + Segment + (1|sample),
                               data = int_data_w_HIM,
                               family = gaussian(), cores = 4,
                               prior = c(prior_b0_depth,
                                         prior_b_depth,
                                         prior_sd_depth,
                                         prior_re_sd_depth),
                               iter = 10000, warmup = 1000,
                               seed = 78)

# plot the posterior distributions of the coefficients
results_depth <- plot_posteriors(int_data_w_HIM, "Segment", depth_model_no_HIM_pred, "Difference to Circle 35 in corrected DPM (cube root)", baseline = "Circle 35", subtract_for_HIM = FALSE, col_by_HIM = TRUE, file_name = "new_version_of_Ines_manuscript/Figure_6.pdf", cube_root = TRUE)
plot_posteriors(int_data_w_HIM, "hostwasp", depth_model_no_HIM_pred, 
                expression("Difference to "~italic("C. s.")~" Kitale -"~italic("B. fusca")~"in corrected DPM (cube root)"), 
                baseline = "B. fusca-C. sesamiae kitale", 
                subtract_for_HIM = FALSE, 
                col_by_HIM = FALSE, 
                file_name = "new_version_of_Ines_manuscript/Figure_7.pdf", 
                cube_root = TRUE,
                pretty_names = c("C. s. Mombasa - B.fusca",
                                 "C. s. Mombasa - S. calamistis",
                                 "C. s. Kitale - S. calamistis",
                                 "C. flavipes - C. partellus",
                                 "C. flavipes - S. frugiperda",
                                 "C. icipe - S. frugiperda",
                                 "C. typhae - S. nonagrioides"))

# check posterior predictions
ppd <- posterior_predict(depth_model_no_HIM_pred)
pdf("new_version_of_Ines_manuscript/Figure_S13A.pdf", width = 4, height = 4)
pp_check(depth_model_no_HIM_pred, ndraws = 1000) + 
  labs(x = "Corrected DPM") + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        text = element_text(family = "sans")) +
  theme(legend.position="none")
dev.off()
pdf("new_version_of_Ines_manuscript/Figure_S13B.pdf", width = 5, height = 4)
plot_ppd_by_segment(ppd, 1000, int_data_w_HIM, "depth_cube_root_cpm",
         "Corrected DPM", logged = FALSE)
dev.off()
######

# play with FASTAs
######
# order segments by depth
ordered <- aggregate(depth_cube_root ~ Segment, data = int_data_w_HIM,
                     FUN = median)
ordered <- ordered[order(ordered[,2], decreasing = TRUE),1]

fastas <- data.frame(wasp = rep(levels(int_data$wasp), each = length(levels(int_data$Segment))),
                     segment = rep(levels(int_data$Segment), length(levels(int_data$wasp))),
                     lengths = rep(0, length(levels(int_data$wasp)) * length(levels(int_data$Segment))),
                     GC = rep(0, length(levels(int_data$wasp)) * length(levels(int_data$Segment))),
                     motif = rep(0, length(levels(int_data$wasp)) * length(levels(int_data$Segment))))
for (w in 1:length(original_wasps)) {
  fasta <- data.frame(readAAStringSet(paste0("kmers/", original_wasps[w], ".fasta"),format="fasta"))
  for (s in rownames(fasta)) {
    curr_length <- nchar(fasta[s,])
    curr_GC <- ((str_count(fasta[s,], "G") + str_count(fasta[s,], "C"))/curr_length) * 100 
    curr_motif <- ((str_count(fasta[s,], "AAAAAAAA"))/curr_length)
    fastas[fastas$wasp == levels(int_data$wasp)[w] & fastas$segment == s, "lengths"] <- curr_length
    fastas[fastas$wasp == levels(int_data$wasp)[w] & fastas$segment == s, "GC"] <- curr_GC
    fastas[fastas$wasp == levels(int_data$wasp)[w] & fastas$segment == s, "motif"] <- curr_motif
  }
}
ordered <- gsub("Segment_", "Circle ", ordered)
fastas$segment <- gsub("Segment_", "Circle ", fastas$segment)
fastas$segment <- factor(fastas$segment, levels = ordered)
pdf("figures/Figure_R8.pdf", width = 7, height = 5)
par(mfrow = c(1,2))
trans_goldenrod <- rgb(244, 196, 48, maxColorValue = 255, alpha = 10)
stripchart(lengths ~ segment, vertical = TRUE,
           method = "jitter", pch = 20, data = fastas, 
           las = 2, cex.axis = 0.7, ylab = "Circle length (nt)")
starting_point <- 1
for (seg in levels(fastas$segment)) {
  if (starting_point%%2 == 0) {
    rect(starting_point - 0.5, 0, starting_point + 0.5, 40000,
         col = alpha("goldenrod1", alpha = 0.4),
         border = NA)
  }
  starting_point <- starting_point + 1
}
stripchart(GC ~ segment, vertical = TRUE,
           method = "jitter", pch = 20, data = fastas, 
           las = 2, ylim = c(28, 37.5), cex.axis = 0.7,
           ylab = "GC (%)")
starting_point <- 1
for (seg in levels(fastas$segment)) {
  if (starting_point%%2 == 0) {
    rect(starting_point - 0.5, 0, starting_point + 0.5, 100,
         col = alpha("goldenrod1", alpha = 0.4),
         border = NA)
  }
  starting_point <- starting_point + 1
}
dev.off()

######
