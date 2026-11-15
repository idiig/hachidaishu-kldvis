# One-shot preprocessing: turns data/reference/kondo2011list.tsv into
# the two lists plot_point_vs_target.R reads directly as files, so
# run_meta_config_plots.R doesn't have to recompute them every run.
# Re-run this (from the visualization root) after editing
# ../data/reference/kondo2011list.tsv:
#   Rscript scripts/visualization/prep_kondo_lists.R
#
# Two outputs:
#   data/reference/kondo2011-annotate-words.txt -- single column, every
#     hl==1 row (expression + stem tokens split on "|" with length>=2),
#     used for plain badges (no type/color).
#   data/reference/kondo2011-highlight.tsv -- term/type two columns,
#     every is_word==1 row, split into two types:
#       "近藤(2011)の抽出語" -- the whole expression
#       "近藤(2011)の抽出語の構成要素" -- stem tokens split on "|" with
#         length>=2 (single-character particle/auxiliary-verb stems are
#         filtered out -- they collide with unrelated entries everywhere)

source("scripts/visualization/text_utils.R")

kondo_expr <- read.table("../data/reference/kondo2011list.tsv", sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8",
                          quote = "", comment.char = "")

hl_rows <- kondo_expr[kondo_expr$hl == 1, ]
annotate_words <- unique(c(trimws(hl_rows$expression), stem_tokens_ge2(hl_rows$stem)))
writeLines(annotate_words, "data/reference/kondo2011-annotate-words.txt", useBytes = TRUE)

iw_rows <- kondo_expr[kondo_expr$is_word == 1, ]
# "構成要素" (component) rows are taken only where stem genuinely
# contains "|" (Kondo split that row into multiple components) -- a row
# where stem equals expression verbatim (no "|", never actually split,
# e.g. "をみなへし" being its own stem) was never decomposed and must
# not also count as a "構成要素", or it would appear under both types
# and (since "構成要素" is rbind'd first, overwritten by whichever type
# plot_point_vs_target.R processes last) get mis-colored.
iw_decomposed <- iw_rows[grepl("|", iw_rows$stem, fixed = TRUE), ]
# "構成要素" listed before "抽出語" -- plot_point_vs_target.R processes
# the highlight table in row order, later rows overwrite earlier ones;
# on the rare chance a component token happens to also be some other
# row's full expression, it should win as the more specific "抽出語"
# classification, not fall back to "構成要素".
highlight <- rbind(
  data.frame(term = stem_tokens_ge2(iw_decomposed$stem),
             type = "近藤(2011)の抽出語の構成要素", stringsAsFactors = FALSE),
  data.frame(term = unique(trimws(iw_rows$expression)),
             type = "近藤(2011)の抽出語", stringsAsFactors = FALSE)
)
write.table(highlight, "data/reference/kondo2011-highlight.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

cat("wrote", length(annotate_words), "annotate words ->",
    "data/reference/kondo2011-annotate-words.txt\n")
cat("wrote", nrow(highlight), "highlight rows ->",
    "data/reference/kondo2011-highlight.tsv\n")
