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

scale_fillcolor_manual <- function(...) {
  return(list(scale_fill_manual(...), scale_color_manual(...)))
}