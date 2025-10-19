#!/usr/bin/env Rscript

###Analyze MIT results from generated sequences to validate expected behavior.

library(ggplot2)
library(MASS)
library(getopt)

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#Use getopt to parse command line options
spec = matrix(c(
  'input', 'i', 1, "character", 
  'output', 'o', 1, "character", 
  'help', 'h', 0, "logical"
), byrow=TRUE, ncol=4)

  opt = getopt(spec)

  if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
  }

##FUNCTIONS##

#This function performs the RCW conservation adjustment on MIT scores
#There are some positions that are more conserved than others, which can lead to
#artificially high MIT scores. This function adjusts for that by normalizing
#the MIT score based on the conservation of the individual positions.
#RCW stands for "Row Column Weighting" - there are other methods to do this as well.
# might be worth exploring in future versions.
rcw_conservation <- function(x) {
  new <- x
  new$RCW_entropy <- rep(0, nrow(x))
  for ( i in 1:nrow(x) ) {
    pos1 <- x[i,1]
    pos2 <- x[i,2]
    
    sum1 <- sum(x$MIT_Score[x$Position_MSA1 == pos1])
    sum2 <- sum(x$MIT_Score[x$Position_MSA2 == pos2])
    #Perform the calculation
    new$RCW_entropy[i] <-  x[i,3]/
    (( (sum1+sum2) - (2*x[i,3]) )/
         ( length(x$MIT_Score[x$Position_MSA1 == pos1]) + length(x$Position_MSA1[x$Position_MSA2 == pos2]) - 2 ))
  }
  return(new)
}

z_score <- function(vect_numbers) {
    # return z-scores vectorized
    return((vect_numbers - mean(vect_numbers, na.rm = TRUE)) / sd(vect_numbers, na.rm = TRUE))
}

main <- function(opt) {
    # Read in the MIT results
    mit_data <- read.csv(opt$input, stringsAsFactors = FALSE)

    mit_data$MIT_Score <- as.numeric(mit_data$MIT_Score)

    #get rid of scores that are with itself
    mit_data <- mit_data[mit_data$Position_MSA1 != mit_data$Position_MSA2, ]

    # Perform RCW conservation adjustment
    mit_data_rcw <- rcw_conservation(mit_data)

    mit_data_rcw <- na.omit(mit_data_rcw)

    # Calculate z-scores
    mit_data_rcw$Z_Score <- z_score(mit_data_rcw$RCW_entropy)

    # Ensure output directory exists
    outdir <- opt$output
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    # Output plots
    png(file.path(outdir, "MIT_distribution.png"))
    print(
        ggplot(mit_data_rcw, aes(x = MIT_Score)) +
        geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
        theme_minimal() +
        labs(title = "MIT Score Distribution", x = "MIT Score", y = "Count")
    )
    dev.off()

    png(file.path(outdir, "RCW_MIT_distribution.png"))
    print(
        ggplot(mit_data_rcw, aes(x = RCW_entropy)) +
        geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
        theme_minimal() +
        labs(title = "MIT Score Distribution RCW corrected", x = "RCW MIT Score", y = "Count")
    )
    dev.off()

    png(file.path(outdir, "MIT_ZScore_heatmap.png"), width = 800, height = 800)
    print(
        ggplot(mit_data_rcw, aes(x = Position_MSA1, y = Position_MSA2, fill = Z_Score)) +
        geom_tile(colour = "white") +
        scale_fill_gradient(high = "firebrick", low = "white") +
        coord_fixed() +
        xlab("Position 1") +
        ylab("Position 2") +
        labs(fill = "zscore(RCW MIT Entropy)") +
        theme_dark()
    )
    dev.off()

    return(0)
}

# Run main only when executed as a script (not when sourced)
if (sys.nframe() == 0L) {
    status <- main(opt)
    quit(save = "no", status = as.integer(status))
}