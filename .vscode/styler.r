if (!("optparse" %in% installed.packages())) install.packages("optparse")
if (!("styler" %in% installed.packages())) install.packages("styler")
library(optparse)
library(styler)

parser <- OptionParser(description = "Style R code")

parser <- add_option(parser,
  "--i",
  action = "store"
)

argv <- parse_args(parser)

style_file(argv$i)
