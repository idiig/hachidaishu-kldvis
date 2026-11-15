# One-shot preprocessing: derives an hl pair list from
# ../data/reference/kondo2011list.tsv's own stem column -- sits alongside
# prep_kondo_lists.R without changing it. Re-run (from the
# visualization root) after editing kondo2011list.tsv:
#   Rscript scripts/visualization/prep_pair_kondo_lists.R
#
# hl definition: the same selection criterion dh-book-draft's
# demo.org:308-313 table used ("whether two entries both belong to the
# same kondo2011 string" is itself the selection rule) -- every row
# whose stem (split on "|") has >= 2 components, all pairwise
# combinations (combn(n, 2)) count as one hl pair. Same decomposition
# logic as dh-book-draft's scripts/hapax_pairs_male.R's covered_keys,
# reusing prep_kondo_lists.R's own nchar>=2 filter (single-character
# particle/auxiliary-verb stems collide everywhere in the entry field;
# this is just used as a badge list, so the cost of noise outweighs the
# cost of is_covered()-style defensive lookups). Not restricted to
# hl==1 -- every row whose stem has >=2 components is included.

source("scripts/visualization/text_utils.R")

kondo_expr <- read.table("../data/reference/kondo2011list.tsv", sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8",
                          quote = "", comment.char = "")

multi <- kondo_expr[grepl("|", kondo_expr$stem, fixed = TRUE), ]

pairs <- do.call(rbind, lapply(multi$stem, function(stem_str) {
  stems <- stem_tokens_ge2(stem_str)
  if (length(stems) < 2) return(NULL)
  combos <- combn(stems, 2, simplify = FALSE)
  do.call(rbind, lapply(combos, function(p) data.frame(w1 = p[1], w2 = p[2])))
}))
pairs <- unique(pairs)
# The "type" column plays the same role as the single-word highlight
# table's own "type" column -- the legend's own label comes from here
# (the data side), not a CLI argument plot_point_vs_target.R would
# need a separate mechanism for.
pairs$type <- "近藤(2011)の抽出語の構成要素対"

write.table(pairs, "data/reference/kondo2011-pair-hl.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

cat(nrow(multi), "multi-component kondo2011 entries ->",
    nrow(pairs), "hl pairs -> data/reference/kondo2011-pair-hl.tsv\n")
