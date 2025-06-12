# R CUSTOM FUNCTIONS / SETTINGS FOR ANALYSES
library(plyr)
library(dplyr)
library(paletteer)

## is.nan function for data frames
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

## notin
`%notin%` <- Negate(`%in%`)

## mode (split ties randomly)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(rank(tabulate(match(x, ux)), ties.method = "random"))]
}

## percentile rank
library(stringr)
percentile_rank <- function(x) 
  ifelse(is.na(x),NA,rank(x)/length(x[is.na(x) == FALSE]))

# function to get closest quantile value
percentile_level <- function(x){
  probs_vec <- quantile(x, probs = seq(0,1,0.001))
  out <- data.frame(orig = x, perc = 0)
  for(i in 1:length(x)){
    abs.diff = abs(x[i] - probs_vec)
    prob.ind = which.min(abs.diff)
    perc.val = as.numeric(str_replace_all(names(probs_vec[prob.ind]), "%", ""))/100
    out[i,2] = perc.val
  }
  return(out[,2])
}

## floor to any interval
floor_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}

# VIF from Zuur (2009)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  # VIF
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Support function for corvif. Will not be called by the user.
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# GRAPHICS
## color blind friendly palette (with grey):
col1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## color_blind_friendly palette (with black): 
col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## larger text size for ggplot
library(ggplot2)
mythemes <- list(
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)))

## phyloseq_to_adonis function from MicrobeMiseq GitHub
phyloseq_to_adonis <- function(physeq, distmat = NULL, dist, formula) {
  
  if(!is.null(distmat)) {
    phydist <- distmat
  } else {
    phydist <- phyloseq::distance(physeq, dist)
  }
  
  metadata <- as(sample_data(physeq), "data.frame")
  
  # Adonis test
  f <- reformulate(formula, response = "phydist")
  adonis.test <- adonis(f, data = metadata)
  print(adonis.test)
  
  # Run homogeneity of dispersion test if there is only 1 variable
  if (grepl("\\+", formula)) {
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test
    )
  } else {
    group <- metadata[ ,formula]
    beta <- betadisper(phydist, group)
    disper.test = permutest(beta)
    print(disper.test)
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test, 
      disper = disper.test
    )
  }
  return (l)
}

# labelling figures in base R
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

# labelling a facetted figure in ggplot2 (https://forum.posit.co/t/how-to-automatically-add-text-annotations-or-tags-outside-of-faceted-plots/13700/3)
tag_facet_flex <- function(p, position = 'both', 
                           tag_fun_top = function(i) letters[i],
                           tag_fun_right = utils::as.roman,
                           x = c(0, 0), y = c(0.5, 1),
                           hjust = c(0, 0), vjust = c(0.5, 1),
                           fontface = c(1, 1), ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout  
  
  if (grepl(position, 'top')) {
    
    lay <- gb$layout$layout
    
    tags_top <- paste0(tag_fun_top(unique(lay$COL)))
    
    tl <- lapply(tags_top, grid::textGrob,
                 x = x[1], y = y[1],
                 hjust = hjust[1], vjust = vjust[1], 
                 gp = grid::gpar(fontface = fontface[1], ...)
    )
    
    g <- ggplot_gtable(gb)
    g <- gtable::gtable_add_rows(g, grid::unit(1, "line"), pos = 0)
    lm <- unique(g$layout[grepl("panel", g$layout$name), "l"])
    g <- gtable::gtable_add_grob(g, grobs = tl, t = 1, l = lm)
    
  } else if (grepl(position, 'right')) {
    
    tags_right <- paste0(open[2], tag_fun_right(unique(lay$ROW)), close[2])
    
    rl <- lapply(tags_right, grid::textGrob,
                 x = x[2], y = y[2],
                 hjust = hjust[2], vjust = vjust[2], 
                 gp = grid::gpar(fontface = fontface[2], ...)
    )
    
    g <- ggplot_gtable(gb)
    l <- unique(g$layout[grepl("panel", g$layout$name), "l"])
    
    wm <- do.call(grid::unit.pmax, lapply(rl, grid::grobWidth))
    g <- gtable::gtable_add_cols(g, wm, pos = max(l))
    t <- unique(g$layout[grepl("panel", g$layout$name), "t"])
    g <- gtable::gtable_add_grob(g, grobs = rl, t = t, l = max(l) + 1)
    g <- gtable::gtable_add_cols(g, unit(2, "mm"), pos = max(l))
    
  } else {
    
    print('Use default tagging option: both top and right sides')
    
    tags_top <- paste0(open[1], tag_fun_top(unique(lay$COL)), close[1])
    tags_right <- paste0(open[2], tag_fun_right(unique(lay$ROW)), close[2])
    
    tl <- lapply(tags_top, grid::textGrob,
                 x = x[1], y = y[1],
                 hjust = hjust[1], vjust = vjust[1], 
                 gp = grid::gpar(fontface = fontface[1], ...)
    )
    
    rl <- lapply(tags_right, grid::textGrob,
                 x = x[2], y = y[2],
                 hjust = hjust[2], vjust = vjust[2], 
                 gp = grid::gpar(fontface = fontface[2], ...)
    )
    
    g <- ggplot_gtable(gb)
    g <- gtable::gtable_add_rows(g, grid::unit(1, "line"), pos = 0)
    l <- unique(g$layout[grepl("panel", g$layout$name), "l"])
    g <- gtable::gtable_add_grob(g, grobs = tl, t = 1, l = l)
    
    wm <- do.call(grid::unit.pmax, lapply(rl, grid::grobWidth))
    g <- gtable::gtable_add_cols(g, wm, pos = max(l))
    t <- unique(g$layout[grepl("panel", g$layout$name), "t"])
    g <- gtable::gtable_add_grob(g, grobs = rl, t = t, l = max(l) + 1)
    g <- gtable::gtable_add_cols(g, unit(2, "mm"), pos = max(l))
  }
  
  if (!is.null(g)) {
    grid::grid.newpage()
    grid::grid.draw(g)    
  }
  
  return(g)
  
}

sample.data.frame <- function(ps) {
  return(as(phyloseq::sample_data(ps), "data.frame"))
}

# vegan cov ellipse (missing from vegan package)
`veganCovEllipse` <- function(cov, center = c(0,0), scale = 1, npoints = 100) {
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o]))
}
