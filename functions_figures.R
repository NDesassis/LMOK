# Auxiliary functions for graphics

#' Initialize the graphic functions for the control of the output (display or files)
#' 
#'@param study.key    A string giving to identify the study
#'@param flag.graphic A boolean value. If TRUE the graphic is computed, else nothing is done.
#'@param flag.file    A boolean value. If TRUE the graphic is store in a file if FALSE the graphic is displayed.
#'@param type.file    A string to specify the file format to store the graphic: png or eps
#'@param dir.fig      A string to specify the directory where the file is written
#'
#'@return A list of functions
#'  - start   to initialize the display/file
#'  - end     to close the file
fig.init <- function(study.key, flag.graphic = TRUE, 
                     flag.file = FALSE, 
                     type.file = "png",
                     dir.fig = "./figures/") {

    fn_start <- function(file_name, radix = "", h = 480, w = 480) {
    if (flag.graphic & flag.file) {
      if(type.file == "png") {
        f_n <- paste0(dir.fig, study.key, "_", radix, "_", file_name, ".png")
        png(filename = f_n, height=h, width=w, units = "px")    
      } else if(type.file == "eps") {
        f_n <- paste0(dir.fig, study.key, "_", radix, "_", file_name, ".eps")
        setEPS()
        postscript(file = f_n)
      }
      print(paste0(" Saving graphic in ", f_n))
    }
    return (flag.graphic)
  }
  
  fn_end <- function() if(flag.graphic & flag.file) dev.off()
  
  list(start = fn_start, end = fn_end)
}

