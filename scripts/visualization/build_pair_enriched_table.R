# 把pair坐标表拼上两侧的freq/dispersion和两侧间的association值，供
# plot_point_vs_target.R的pair模式使用。只做join，不重新计算KLD/计数。
#
# 用法(key=value风格):
#   Rscript build_pair_enriched_table.R \
#     pair_input=data/processed/point-vs-target/pair-entry-vs-gender-kokin.tsv \
#     single_input_w1=data/processed/point-vs-target/entry-vs-gender-kokin.tsv \
#     single_input_w2=data/processed/point-vs-target/entry-vs-gender-kokin.tsv \
#     assoc_input=data/processed/point-vs-target/pair-entry-entry-association-kokin.tsv \
#     output=data/processed/point-vs-target/pair-entry-vs-gender-kokin-enriched.tsv

parse_cli_flags <- function(args) {
  kv <- strsplit(args, "=", fixed = TRUE)
  setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))
}

flag <- function(flags, name, default = NULL) {
  if (name %in% names(flags)) flags[[name]] else default
}

read_tsv <- function(path) {
  read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             quote = "", fileEncoding = "UTF-8", na.strings = c("NA", ""))
}

# visprepare.py现在把term/doc两种frequency_type合并进同一个文件(point
# 重复一份per frequency_type)——total_freq/n_target/n_not_target/assoc_*
# 这几列的数值真的随frequency_type变(dispersion虽然不变但统一走同一套
# 查表逻辑，不单独开一条路径)，所以查表key必须带上frequency_type，不能
# 只按point名字查(否则多份重复point行里随机命中哪一份)
point_lookup <- function(table, value_col) {
  key <- paste(table$point, table$frequency_type, sep = "\t")
  setNames(table[[value_col]], key)
}

# keyness(单词自身针对某个target的signed KLD)跟freq1/dispersion1不一样，
# 真的随kld_variant变(obs/nofreq/norm/keyassoc是4种不同的计算)——不能
# 只按(point, frequency_type)查，会在visprepare.py把多个variant合并进
# 同一份单词表之后随机命中其中一个variant的行(2026-11-07新增，给
# filter_by_node_keyness/draft.org:1211-1212用)
point_lookup_by_variant <- function(table, value_col) {
  # Category tables (e.g. visprepare's category_point_output) carry no
  # kld_variant/keyness: every lookup is then NA rather than an error.
  if (is.null(table$kld_variant) || is.null(table[[value_col]])) {
    return(setNames(numeric(0), character(0)))
  }
  key <- paste(table$point, table$frequency_type, table$kld_variant, sep = "\t")
  setNames(table[[value_col]], key)
}

load_single_word_tables <- function(flags) {
  # entry×entry这种field1==field2的pair两侧共用同一份表，指定同一个路径时
  # 只读一次，避免重复I/O
  w1_path <- flag(flags, "single_input_w1", flag(flags, "single_input"))
  w2_path <- flag(flags, "single_input_w2", flag(flags, "single_input"))
  w1 <- read_tsv(w1_path)
  w2 <- if (identical(w2_path, w1_path)) w1 else read_tsv(w2_path)
  list(w1 = w1, w2 = w2)
}

report_missing_join <- function(d, cols, label) {
  missing <- Reduce(`|`, lapply(cols, function(col) is.na(d[[col]])))
  cat(sum(missing), "/", nrow(d), label, "\n")
}

join_single_word_stats <- function(d, single) {
  w1_key <- paste(d$w1, d$frequency_type, sep = "\t")
  w2_key <- paste(d$w2, d$frequency_type, sep = "\t")
  d$freq1 <- point_lookup(single$w1, "total_freq")[w1_key]
  d$freq2 <- point_lookup(single$w2, "total_freq")[w2_key]
  report_missing_join(d, c("freq1", "freq2"), "pairs missing single-word freq lookup")
  d$dispersion1 <- point_lookup(single$w1, "dispersion")[w1_key]
  d$dispersion2 <- point_lookup(single$w2, "dispersion")[w2_key]
  report_missing_join(d, c("dispersion1", "dispersion2"), "pairs missing single-word dispersion lookup")
  # keyness1/keyness2(2026-11-07新增)：pair片侧(w1/w2)各自的单词单位
  # keyness，供plot_point_vs_target.R的filter_by_node_keyness用(draft.
  # org:1211-1212「それぞれの語の女性歌への連関が0.8以下のものに限定する」)。
  # 跟freq1/dispersion1不同，查表key必须带kld_variant(见point_lookup_
  # by_variant自己的注释)
  w1_variant_key <- paste(d$w1, d$frequency_type, d$kld_variant, sep = "\t")
  w2_variant_key <- paste(d$w2, d$frequency_type, d$kld_variant, sep = "\t")
  d$keyness1 <- point_lookup_by_variant(single$w1, "keyness")[w1_variant_key]
  d$keyness2 <- point_lookup_by_variant(single$w2, "keyness")[w2_variant_key]
  report_missing_join(d, c("keyness1", "keyness2"), "pairs missing single-word keyness lookup")
  d
}

resolve_association <- function(d, assoc) {
  # association表每个无序pair、每个frequency_type各存一份、w1/w2固定顺序，
  # 查询方向可能相反，查不到就反过来查一次、对调assoc_w1_to_w2/
  # assoc_w2_to_w1——key必须带frequency_type，否则term/doc两份数值不同的
  # 行会随机互相覆盖
  key <- paste(assoc$w1, assoc$w2, assoc$frequency_type, sep = "\t")
  a2b <- setNames(assoc$assoc_w1_to_w2, key)
  b2a <- setNames(assoc$assoc_w2_to_w1, key)
  fwd_key <- paste(d$w1, d$w2, d$frequency_type, sep = "\t")
  rev_key <- paste(d$w2, d$w1, d$frequency_type, sep = "\t")
  v12 <- a2b[fwd_key]; v21 <- b2a[fwd_key]
  needs_reverse <- is.na(v12) | is.na(v21)
  v12[needs_reverse] <- b2a[rev_key[needs_reverse]]
  v21[needs_reverse] <- a2b[rev_key[needs_reverse]]
  list(w1_to_w2 = unname(v12), w2_to_w1 = unname(v21))
}

join_association <- function(d, assoc_path) {
  resolved <- resolve_association(d, read_tsv(assoc_path))
  d$assoc_w1_to_w2 <- resolved$w1_to_w2
  d$assoc_w2_to_w1 <- resolved$w2_to_w1
  report_missing_join(d, c("assoc_w1_to_w2", "assoc_w2_to_w1"), "pairs missing association data")
  d
}

write_enriched_table <- function(d, output_path) {
  write.table(d, output_path, sep = "\t", row.names = FALSE, quote = FALSE,
              fileEncoding = "UTF-8", na = "")
  cat("wrote", nrow(d), "rows,", ncol(d), "cols ->", output_path, "\n")
}

flags <- parse_cli_flags(commandArgs(trailingOnly = TRUE))
single <- load_single_word_tables(flags)
# w1/w2直接用pair_input自己的原始列(extract_point_vs_target.py按pair内部"|"
# 分隔符拆出来的)，不从显示名重新切分——按"・"反拆会被entry自身内部的"・"切错
d <- read_tsv(flag(flags, "pair_input"))
d <- join_single_word_stats(d, single)
d <- join_association(d, flag(flags, "assoc_input"))
write_enriched_table(d, flag(flags, "output"))
