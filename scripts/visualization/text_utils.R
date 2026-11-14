stem_tokens_ge2 <- function(stem_col) {
  stems <- trimws(unlist(strsplit(stem_col, "\\|")))
  unique(stems[nchar(stems) >= 2])
}
