#' bundle up some ugly plotting utility code
#' 
#' @param p, a ggplot / grob object
#' @param heightdefault,
#' @param widthdefault, numeric; the height & width defaults (for when shell environment WIDTH & HEIGHT are unset)
#' @param args, the command line args, with the last item being the target file
#' @param units, the units for h&w (both for the defaults or the shell env. values)
#' @param dpi, the resolution (dots per inch)
#' @return same as \code{dev.off()}
#' @example 
#' args <- commandArgs(trailingOnly = TRUE)
#' p <- ggplot(...) // some ggplot
#' plotutil(p, 5, 10, args)
#' 
#' @export
plotutil <- function(p, heightdefault, widthdefault, args="Rplot.png", units="in", res=300) {
  targetfile <- tail(args,1)
  ispng <- grepl("png$", targetfile)
  istiff <- !ispng && grepl("tiff$", targetfile)
  ispdf <- grepl("pdf$", targetfile)
  figdim <- within(as.list(Sys.getenv(c("WIDTH","HEIGHT"), unset = NA)),{
    if (is.na(HEIGHT)) HEIGHT <- heightdefault
    if (is.na(WIDTH)) WIDTH <- widthdefault
    HEIGHT <- as.numeric(HEIGHT)
    WIDTH <- as.numeric(WIDTH)
  })
  with(figdim,{
    if (ispng) {
      png(targetfile, width = WIDTH, height = HEIGHT, units = units, res = res)
    } else if (istiff) {
      tiff(targetfile, width = WIDTH, height = HEIGHT, units = units, res = res, compression = "lzw+p", type = "cairo")
    } else if (ispdf) {
      pdf(targetfile, width = WIDTH, height = HEIGHT)
    } else stop("did not understand device...")
    print(p)
    dev.off()
  })
}

dbutil <- function(dbfile, sql) {
  drv = dbDriver("SQLite")
  db = dbConnect(drv, dbfile, flags=SQLITE_RO)
  res <- data.table(dbGetQuery(db, sql))
  dbDisconnect(db)
  return(res)
}

#' Convenience constructor for matched color + fill ggplot::scale
#'
#' @param ..., any arguments to color + fill manual scales
#' @return a list of ggplot scale objects
#' OBE due to aesthetics arg in scale_fill|color?
scale_fillcolor_manual <- function(...) list(
  ggplot::scale_fill_manual(...),
  ggplot::scale_color_manual(...)
)

#' Convenience method providing different y-axis labels on facet plots
#' with different y scales (using facet labels as the scale labels)
#'
#' MUST BE USED AFTER ANY COMPLETE THEME (e.g., `theme_minimal()`), I.E.
#' `gg` OBJECTS THAT RESET THEME ELEMENTS
#'
#' @param ..., typical arguments to facet_grid *except* `scales` and `switch`
facet_grid_freey <- function(...) list(
  facet_grid(..., scales = "free_y", switch = "y"),
  theme(axis.title.y = element_blank(), strip.placement = "outside")
)

geom_limits <- function(dt) geom_blank(
  mapping = aes(
    shape=NULL, fill=NULL, color=NULL,
    linetype=NULL, size=NULL, group=NULL, alpha=NULL
  ),
  data = dt
)

#' Convenience method for calculating quantiles in data.table context
#'
#' @param probs, same as \code{quantile} `probs` argument (though can be a list).
#' @param val, the data.table column to be \code{quantile}'d
#' @param pnames, resulting column names in data.table; defaults to names of `probs`
#' @param na.rm, same as \code{quantile}; defaults to `TRUE`
#' @return same as \code{quantile}, but as a list with appropriate names
#' @example 
#' dt <- data.table(x=runif(100), y=sample(2, size=100, repl=T))
#' dt[, dtquantiles(.(med=.5), x), by=y]
#' dt[, dtquantiles(.(med=.5), x)]
#' # compare to
#' median(dt$x)
dtquantiles <- function(probs, val, pnames = names(probs), na.rm = T, mn = F) {
  qs <- quantile(val, probs = unlist(probs), na.rm = na.rm)
  if (mn) {
    qs <- c(qs, mean(val, na.rm = na.rm))
    pnames <- c(pnames, "mean")
  }
  names(qs) <- pnames
  as.list(qs)
}

#' Factory-pattern function for making reuseable (gg|cow)plot scale generator
#' 
#' Returns a function that can be invoked to make `scale_%_manual` elements with several
#' defaults set.
#'
#' @param scale, what aesthetic being scaled; defaults to "color", def arg list shows
#'  which scales can be generated
#' @param defname, simple string (or any other legal `name` arg to `scale_manual``);
#'  default name for scale
#' @param deflabels, named vector of characters (or any other legal `labels` arg to `scale_manual``);
#'  default labels
#' @param defvalues, named vector (or any other legal `values` arg to `scale_manual``); default values
#' @return function(name=defname, labels=deflabels, values=defvalues, ...) that can be invoked to produce
#'  scale_%_manual
#' @example 
#' require(ggplot2)
#' require(data.table)
#' scale_color_things <- scale_generator("color", "Things",
#'   c(thing1="One", thing2="Two"), c(thing1="red",thing2="blue")
#' )
#' p.dt <- data.table(x=1:10, y=1:10, thing=sample(c("thing1","thing2"), 10, rep=T))
#' p.dt[thing == "thing1"]
#' p.dt[thing == "thing2"]
#' ggplot(p.dt) +
#'  aes(x=x,y=y,color=thing) + geom_point() +
#'  scale_color_things()
scale_generator <- function(
  scale = c("color", "fill", "linetype", "size", "shape", "alpha")[1],
  defname, deflabels, defvalues
) {
  scalef <- switch(scale,
    color = ggplot2::scale_color_manual,
    fill = ggplot2::scale_fill_manual,
    linetype = ggplot2::scale_linetype_manual,
    size = ggplot2::scale_size_manual,
    shape = ggplot2::scale_shape_manual,
    alpha = ggplot2::scale_alpha_manual
  )
  return(function(
    name=defname, labels=deflabels, values=defvalues,
    ...
  ) {
    otherargs <- list(...)
    if(any(grepl("breaks",names(otherargs)))) labels <- labels[as.character(otherargs$breaks)]
    scalef(
      name = name, labels = labels,
      values = values, ...
    )
  })
}

#' Convenience method for assembling tabularx table innards
#'
#' @param dt, a data.table with rows already generated (via, e.g., `sprintf(...)`)
#' @param whchcol, which column has the rows (defaults to "rows", but can also be a number)
#' @return a string suitable for `cat` to a file
tabularxinnards <- function(dt, whchcol = "rows") {
  subdt <- dt[,.SD,.SDcols=whchcol]
  firstrow <- paste0(subdt[1],"\\T ")
  lastrow <- paste0(subdt[.N],"\\B ")
  return(paste0(sprintf("%s\\\\\n",
    c(firstrow, unlist(subdt[-c(1,.N)]), lastrow)
  ), collapse = ""))
}

logistic <- function(x, x0=0, k=1, L=1) L/(1+exp(-k*(x-x0)))
