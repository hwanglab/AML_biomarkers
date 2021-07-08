WriteInvocation <- function(argv, output_path) {
  argv2 <- argv[2:length(argv)]

  header <- commandArgs()[6]

  args_used <- paste(names(argv2), argv2, sep = ": ")
  header2 <- paste(header, format(Sys.time()), sep = ": ")
  header2 <- stringr::str_remove(header2, "--file=")
  header3 <- paste0(header2, "------------------")

  cat(
    header3, args_used, "\n",
    file = here(output_path),
    sep = "\n",
    append = TRUE
  )
}
