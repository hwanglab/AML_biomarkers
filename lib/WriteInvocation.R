WriteInvocation <- function(argv, output_path) {
  argv2 <- argv

  header <- commandArgs() %>%
    str_subset("radi") %>%
    #stringr::str_remove("--file=") %>%
    basename() #%>%
    #str_remove(".R")

  args_used <- paste(names(argv2), argv2, sep = ": ")

  cat(
    args_used, "\n",
    file = here(output_path, header),
    sep = "\n"
  )
}

TestInvocation <- function(argv, output_path) {
  argv2 <- argv[2:length(argv)]

  header <- args %>% str_subset("--file=")

  header <- header %>%
    str_subset("--file=") %>%
    stringr::str_remove("--file=") %>%
    basename() %>%
    str_remove(".R")
  
  exarg <- read_delim(here(output_path, header), delim = ": ", col_names = c("arg", val))
  exarg$val[is.na(exarg$val)] <- ""
  argt <- as.list(exarg$val)
  names(argt) <- exarg$arg

  test_res <- all.equal(argv2, argt)

  return(isTRUE(is.logical(test_res)))
}