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
  mapping = aes(color=NULL, linetype=NULL, size=NULL, group=NULL, alpha=NULL),
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
dtquantiles <- function(probs, val, pnames = names(probs), na.rm = T) {
  qs <- quantile(val, probs = unlist(probs), na.rm = na.rm)
  names(qs) <- pnames
  as.list(qs)
}