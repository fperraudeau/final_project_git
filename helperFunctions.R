RotateImageVector <- function(image.vec) {
  # Rotate a raw image vector into the
  # correct orientation for viewing.
  
  # There are this many rows and columns,
  # and transpose and re-ordering gives it the correct
  # orientation.
  kRows <- 128
  return(t(matrix(image.vec, nrow = kRows)[kRows:1,]))
}


ReadImage <- function(number) {
  # Read a single line from fit_stim.csv into
  # a matrix that can be plotted as an image.
  #
  # Args:
  #  - number: Which image number to load, starting
  #            with 1.
  #
  # Returns:
  #  - A matrix containing the raw image pixels.
  
  # The first line contains column names, so we always
  # skip at least one line.
  img <- scan("fit_stim.csv",
              skip = number, nlines = 1, sep=",")
  
  return(RotateImageVector(img))
  
}