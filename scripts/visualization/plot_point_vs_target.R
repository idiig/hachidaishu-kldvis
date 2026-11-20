
library(ggrepel)


.compute_axis_ranges <- function(d, nofreq) {
  if (nofreq) {
    list(xr = c(0, 1), yr = c(0, 1), xpad = 0.03, ypad = 0.03, neutral = 0.5)
  } else {
    list(xr = c(-1, 1), yr = c(0, 1), xpad = 0.06, ypad = 0.03, neutral = 0)
  }
}

.axis_break_transform <- function(data_max, axis_max = 1, break_frac = 0.5, main_frac = 0.82) {
  if (!is.finite(data_max) || data_max <= 0) return(NULL)
  gap <- axis_max - data_max
  if (gap / axis_max < break_frac) return(NULL)
  main_end <- min(data_max * 1.15, axis_max * 0.999)
  fwd <- function(x) {
    ifelse(x <= main_end,
           x / main_end * main_frac,
           main_frac + (x - main_end) / (axis_max - main_end) * (1 - main_frac))
  }
  list(fwd = fwd, main_end = main_end, main_frac = main_frac, axis_max = axis_max)
}

.apply_nofreq_breaks <- function(d, nofreq) {
  if (!nofreq) return(list(d = d, x_break = NULL, y_break = NULL))
  x_break <- .axis_break_transform(max(d$frequency))
  y_break <- .axis_break_transform(max(d$keyness))
  if (!is.null(x_break)) d$frequency <- x_break$fwd(d$frequency)
  if (!is.null(y_break)) d$keyness <- y_break$fwd(d$keyness)
  list(d = d, x_break = x_break, y_break = y_break)
}

.draw_axis_break_mark <- function(pos_frac, side, usr) {
  if (side == 1) {
    xm <- usr[1] + pos_frac * (usr[2] - usr[1])
    y0 <- usr[3]
    dy <- (usr[4] - usr[3]) * 0.018
    dx <- (usr[2] - usr[1]) * 0.008
    segments(xm - dx, y0 - dy, xm + dx, y0 + dy, xpd = NA, lwd = 1.5)
    segments(xm + dx, y0 - dy, xm + 3 * dx, y0 + dy, xpd = NA, lwd = 1.5)
  } else {
    ym <- usr[3] + pos_frac * (usr[4] - usr[3])
    x0 <- usr[1]
    dx <- (usr[2] - usr[1]) * 0.018
    dy <- (usr[4] - usr[3]) * 0.008
    segments(x0 - dx, ym - dy, x0 + dx, ym + dy, xpd = NA, lwd = 1.5)
    segments(x0 - dx, ym + dy, x0 + dx, ym + 3 * dy, xpd = NA, lwd = 1.5)
  }
}

.axis_ticks_with_break <- function(brk) {
  ticks <- pretty(c(0, brk$main_end))
  ticks <- ticks[ticks <= brk$main_end]
  ticks <- unique(c(ticks, brk$axis_max))
  list(at = brk$fwd(ticks), labels = format(ticks, trim = TRUE))
}

.open_canvas <- function(output, canvas_size, width_multiplier = 1, height_multiplier = 1) {
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  canvas_scale <- canvas_size / 300
  panel_h_in <- canvas_size * height_multiplier / 72
  width_in <- canvas_size * width_multiplier / 72
  # PNG output -- a raster device's file size/render cost scales with
  # PIXEL count, not point count, unlike svg()'s one vector element
  # per point. A pair dataset's ~300k-point scatter makes svg() both
  # slow to write and (once embedded inline, see the Shiny app's own
  # output$preview) slow for a browser to parse -- confirmed 2026-10-17,
  # a 300k-point render was a 100MB+ SVG. This project's own LaTeX
  # pipeline only ever consumes a real .svg (see dh-book-draft/
  # CLAUDE.md), so this is purely an additional OUTPUT format this
  # function now also accepts -- any caller that still passes a .svg
  # OUTPUT gets the exact same vector file as before.
  if (grepl("\\.png$", output, ignore.case = TRUE)) {
    png(output, width = width_in, height = panel_h_in, units = "in",
        res = 150, type = "cairo", family = getOption("kld.plot.font_family", "Hiragino Mincho ProN"))
  } else {
    svg(output, width = width_in, height = panel_h_in,
        family = getOption("kld.plot.font_family", "Hiragino Mincho ProN"))
  }
  # Bottom/left margins grew from 3.2/3.5 (2026-10-16) to fit the two
  # new rows .draw_axis_direction_arrows() draws between the tick
  # labels and the (now correspondingly pushed-out, see .draw_axes's
  # own title_line) axis title -- a drawn double-arrow row, then an
  # end-text-label row below it. Both ends of this change read the
  # same canvas_scale, so the two stay proportioned together at any
  # canvas_size.
  par(mar = c(5.0, 5.0, 1, 1) * canvas_scale, mgp = c(2, 0.6, 0) * canvas_scale)
  par(cex.axis = canvas_scale, cex.lab = canvas_scale)
  list(canvas_scale = canvas_scale, panel_h_in = panel_h_in)
}

.THRESHOLD_OPS <- c(gte = ">=", lte = "<=", lt = "<")

.resolve_threshold_op <- function(op) {
  if (is.null(op) || is.na(op)) return("gte")
  norm <- tolower(trimws(op))
  norm <- switch(norm, ">=" = "gte", "<=" = "lte", "<" = "lt",
                 "ge" = "gte", "le" = "lte", norm)
  if (!norm %in% names(.THRESHOLD_OPS)) {
    stop("threshold op must be one of gte/lte/lt (or >=/<=/<, ge/le), got: ", op)
  }
  norm
}

.threshold_mask <- function(d, cols, threshold, op = "gte", use_abs = FALSE) {
  if (is.null(threshold) || !all(vapply(cols, function(c) !is.null(d[[c]]), logical(1)))) {
    return(rep(TRUE, nrow(d)))
  }
  op <- .resolve_threshold_op(op)
  Reduce(`&`, lapply(cols, function(c) {
    v <- d[[c]]
    if (use_abs) v <- abs(v)
    if (identical(op, "gte")) !is.na(v) & v >= threshold else if (identical(op, "lte")) {
      !is.na(v) & v <= threshold
    } else {
      !is.na(v) & v < threshold
    }
  }))
}

.filter_by_columns <- function(d, cols, threshold, op = "gte", use_abs = FALSE) {
  d[.threshold_mask(d, cols, threshold, op, use_abs), ]
}

.SELECTION_FILTER_COLS <- list(
  edge = c("assoc_w1_to_w2", "assoc_w2_to_w1"),
  dispersion = "dispersion",
  node_dispersion = c("dispersion1", "dispersion2"),
  frequency = "total_freq",
  node_frequency = c("freq1", "freq2"),
  node_keyness = c("keyness1", "keyness2"),
  node_signed_keyness = c("keyness1", "keyness2")
)

# どのfilterがabs(v)を取ってから比較するか——今のところnode_keynessだけ
# (符号付きKLDの「絶対値」を見る、draft.org:1211-1212)。ここに無いkeyは
# FALSE扱い(.threshold_mask自身のデフォルトと同じ)。
.SELECTION_FILTER_ABS <- c(node_keyness = TRUE)

.selection_filter_spec <- function(edge = NULL, edge_op = "gte",
                                    dispersion = NULL, dispersion_op = "gte",
                                    node_dispersion = NULL, node_dispersion_op = "gte",
                                    frequency = NULL, frequency_op = "gte",
                                    node_frequency = NULL, node_frequency_op = "gte",
                                    node_keyness = NULL, node_keyness_op = "gte",
                                    node_signed_keyness = NULL,
                                    node_signed_keyness_op = "gte") {
  list(
    edge = list(value = edge, op = edge_op),
    dispersion = list(value = dispersion, op = dispersion_op),
    node_dispersion = list(value = node_dispersion, op = node_dispersion_op),
    frequency = list(value = frequency, op = frequency_op),
    node_frequency = list(value = node_frequency, op = node_frequency_op),
    node_keyness = list(value = node_keyness, op = node_keyness_op),
    node_signed_keyness = list(value = node_signed_keyness, op = node_signed_keyness_op)
  )
}

.selection_filter_idx <- function(d, spec) {
  mask <- rep(TRUE, nrow(d))
  for (key in names(.SELECTION_FILTER_COLS)) {
    f <- spec[[key]]
    if (is.null(f) || is.null(f$value)) next
    use_abs <- key %in% names(.SELECTION_FILTER_ABS) && isTRUE(.SELECTION_FILTER_ABS[[key]])
    mask <- mask & .threshold_mask(d, .SELECTION_FILTER_COLS[[key]], f$value, f$op,
                                    use_abs = use_abs)
  }
  which(mask)
}

.FILTER_DESC_LABEL <- c(edge = "連関", dispersion = "分散度", node_dispersion = "語の分散度",
                         frequency = "頻度", node_frequency = "語の頻度",
                         node_keyness = "語自身の|KLD|",
                         node_signed_keyness = "語自身のKLD")
.FILTER_DESC_OP_SYMBOL <- c(gte = "≥", lte = "≤", lt = "<")

.inset_filter_desc <- function(spec) {
  parts <- character(0)
  for (key in names(.FILTER_DESC_LABEL)) {
    f <- spec[[key]]
    if (!is.null(f) && !is.null(f$value)) {
      op_sym <- .FILTER_DESC_OP_SYMBOL[[.resolve_threshold_op(f$op)]]
      parts <- c(parts, sprintf("%s %s %s", .FILTER_DESC_LABEL[[key]], op_sym,
                                 format(f$value, trim = TRUE)))
    }
  }
  if (length(parts) == 0) return("(フィルタなし)")
  paste(parts, collapse = "、")
}

filter_by_edge_strength <- function(d, threshold, op = "gte") {
  .filter_by_columns(d, c("assoc_w1_to_w2", "assoc_w2_to_w1"), threshold, op)
}

filter_by_dispersion <- function(d, threshold, op = "gte") {
  .filter_by_columns(d, "dispersion", threshold, op)
}

filter_by_node_dispersion <- function(d, threshold, op = "gte") {
  .filter_by_columns(d, c("dispersion1", "dispersion2"), threshold, op)
}

filter_by_frequency <- function(d, threshold, op = "gte") {
  # dispersionとは別の判断基準——dispersionは「分布が1-2首に偏っている」
  # を捉える、これはそもそも「出現回数が少なすぎて信頼できない」を捉える
  # (2026-10-15、ユーザーと確認済み：両方残す、片方をもう片方に差し替える
  # のではない)。word/pair両版表ともtotal_freq列を持つのでpair限定の
  # ガードは不要
  .filter_by_columns(d, "total_freq", threshold, op)
}

filter_by_node_frequency <- function(d, threshold, op = "gte") {
  # pair only, mirrors filter_by_node_dispersion: filter_by_frequency
  # looks at the pair's combined total_freq, this looks at each side
  # entry's own freq1/freq2 -- a frequent pair can still have one rare
  # side word. Both sides must satisfy the threshold. word table has no
  # freq1/freq2, no-op
  .filter_by_columns(d, c("freq1", "freq2"), threshold, op)
}

filter_by_node_keyness <- function(d, threshold, op = "gte") {
  # pair only, mirrors reproduction/pair_feature_list_reproduce.org's
  # own filter_by_node_keyness (draft.org:1211-1212): a word-pair whose
  # OWN single-word KLD toward the target already exceeds the threshold
  # doesn't contribute new information beyond what the single-word
  # analysis already captured, regardless of the pair's combined
  # keyness -- so both sides' own keyness (signed, hence use_abs=TRUE)
  # must stay within the threshold. word table has no keyness1/
  # keyness2, no-op.
  .filter_by_columns(d, c("keyness1", "keyness2"), threshold, op, use_abs = TRUE)
}

# 図中の小さい文字(badge・mean/median・方向矢印の文字)をまとめて
# 拡大する倍率(2026-11-14、draft-short.org用)。draw_point_vs_target()の
# text_scale引数がoptionsに入れる。既定1=従来どおり。凡例は位置を
# インチ固定で並べているため対象外(拡大すると文字同士が重なる)。
.text_scale <- function() getOption("kld.plot.text_scale", 1)

.draw_axes <- function(d, ranges, nofreq, axis_label, point_label, x_break = NULL, y_break = NULL,
                        axis_col = "black", kld_variant = NULL, label_lang = "ja",
                        guides = TRUE) {
  xr <- ranges$xr; yr <- ranges$yr
  plot.new()
  plot.window(xlim = c(xr[1] - ranges$xpad, xr[2] + ranges$xpad),
              ylim = c(yr[1] - ranges$ypad, yr[2] + ranges$ypad))

  kld_sub <- if (nofreq) "nofreq" else "signed"
  # Read-only (never par(mgp = ...) here) -- .draw_axes() runs twice
  # per device in inset mode (source panel then zoomed panel, see
  # .render_with_inset), sharing one .open_canvas() call; mutating
  # par("mgp") here would compound on the second call instead of
  # giving both panels the same offset. title_line pushes the axis
  # title out past where .draw_axis_direction_arrows() now puts its own
  # arrow AND that arrow's own end-text-label row (between the tick
  # labels and the title) -- both read the same base_mgp so they
  # always stay in the same relative order regardless of canvas_scale.
  base_mgp <- par("mgp")
  title_line <- base_mgp[1] + 1.4
  if (nofreq) {
    if (!is.null(x_break)) {
      xt <- .axis_ticks_with_break(x_break)
      axis(1, at = xt$at, labels = xt$labels, col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    } else {
      axis(1, at = pretty(xr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    if (!is.null(y_break)) {
      yt <- .axis_ticks_with_break(y_break)
      axis(2, at = yt$at, labels = yt$labels, col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    } else {
      axis(2, at = pretty(yr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"),
          ylab = bquote(plain(.(.ui_label("dispersion", label_lang))) * "(" * .(point_label) * ")"),
          line = title_line)
    usr <- par("usr")
    if (!is.null(x_break)) .draw_axis_break_mark(x_break$main_frac, side = 1, usr = usr)
    if (!is.null(y_break)) .draw_axis_break_mark(y_break$main_frac, side = 2, usr = usr)
  } else {
    if (!isTRUE(getOption("kld.plot.hide_x_ticks", FALSE))) {
      axis(1, at = pretty(xr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    if (!isTRUE(getOption("kld.plot.hide_y_ticks", FALSE))) {
      axis(2, at = pretty(yr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    if (!isTRUE(getOption("kld.plot.hide_x_title", FALSE))) {
      title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"),
            line = title_line)
    }
    if (!isTRUE(getOption("kld.plot.hide_y_title", FALSE))) {
      title(ylab = bquote(plain(.(.ui_label("frequency", label_lang))) * "(" * .(point_label) * ")"),
            line = title_line)
    }
  }
  # guides = FALSE (zoomed inset panels) skips the neutral/median/mean
  # reference lines and their labels.
  if (!guides) return(invisible())
  neutral_x <- if (!is.null(x_break)) x_break$fwd(ranges$neutral) else ranges$neutral
  neutral_y <- if (!is.null(y_break)) y_break$fwd(ranges$neutral) else ranges$neutral
  abline(h = neutral_y, v = neutral_x, col = "#CCCCCC", lty = 2)
  y_median <- median(d$keyness)
  y_mean <- mean(d$keyness)
  abline(h = y_median, v = median(d$frequency), col = "#AAAAAA", lty = 1)
  # Mean line + text labels for both (2026-10-15, user-requested) --
  # the Y axis's own median/mean only (the X-axis median line above
  # stays unlabeled, out of scope here). Placed at xr[2] (the data
  # range's own right edge, clipped to the panel, xpd=FALSE) rather
  # than the canvas edge -- the color legend/size legend already own
  # the margin past the panel (see .draw_legend/.draw_size_legend),
  # text out there would collide with them.
  abline(h = y_mean, col = "#888888", lty = 3)
  text(xr[2], y_median, labels = sprintf("median %.3f", y_median),
       col = "#AAAAAA", cex = 0.55 * .text_scale(), adj = c(1, -0.3))
  text(xr[2], y_mean, labels = sprintf("mean %.3f", y_mean),
       col = "#888888", cex = 0.55 * .text_scale(), adj = c(1, -0.3))
}

# 点の大きさは2026-11-05以降dispersion基準(draw_point_vs_target()の冒頭で
# d$node_size <- d$dispersionと差し替え済み、.draw_size_legend自身の
# コメント参照)——以前はn_target(頻度)のlog2正規化だったが、dispersionは
# もともと表生成時に0〜1に近い範囲の指標として計算済みの列なので、Rの
# 描画側でさらに正規化し直す必要はなく、そのままsqrt(面積比例)の式に
# 差し込むだけでよい。0.7/0.3/1.1は2026-09-26確認済みの大きさ係数、
# 変更なし(入力の意味が変わっただけで、面積比例の式自体はそのまま)。
.size_scale <- function(node_size) 0.7 * (0.3 + 1.1 * sqrt(node_size))

.capped_zoom_factor <- function(zoom_factor, cap = 4) min(zoom_factor, cap)

field_shape_by_type <- function() c(entry = 1, concept = 0, topic = 5)

field_label_by_type <- function() c(entry = "語", concept = "概念", topic = "部立")


.abs_ramp <- function(vals, abs_range, light_col, dark_col) {
  norm <- if (diff(abs_range) > 0) {
    (abs(vals) - abs_range[1]) / diff(abs_range)
  } else {
    rep(0.5, length(vals))
  }
  ramp <- colorRampPalette(c(light_col, dark_col))(101)
  ramp[pmin(100, pmax(1, round(norm * 100))) + 1]
}

.linear_ramp <- function(vals, val_range, lo_col, hi_col) {
  norm <- if (diff(val_range) > 0) (vals - val_range[1]) / diff(val_range) else rep(0.5, length(vals))
  ramp <- colorRampPalette(c(lo_col, hi_col))(101)
  ramp[pmin(100, pmax(0, round(norm * 100))) + 1]
}

category_color_palettes <- function() list(
  c(light = "#F4B6AC", dark = "#7A0C00"),
  c(light = "#C9E9CE", dark = "#164A22"),
  c(light = "#E5D4F0", dark = "#4A1B6B"),
  c(light = "#FCE8B2", dark = "#8A5A00"),
  c(light = "#B8E2E8", dark = "#0D4F5C")
)
default_color_palette <- function() c(light = "#E5E5E5", dark = "#9E9E9E")
neutral_color_palette <- function() c(light = "#D9D9D9", dark = "#1A1A1A")
gradient_color_palette <- function() c(lo = "#E5E5E5", hi = "#1A1A1A")

.compute_point_colors <- function(d, highlight, highlight_pairs = NULL,
                                   annotate_point_exclude = NULL,
                                   highlight_frontier = FALSE,
                                   frontier_n = 15,
                                   frontier_candidate_idx = NULL,
                                   category_palettes = category_color_palettes(),
                                   default_palette = default_color_palette(),
                                   neutral_palette = neutral_color_palette(),
                                   gradient_palette = gradient_color_palette(),
                                   unselected_alpha = 0.05) {
  # Shade now tracks dispersion_difference (D_T - D_R, signed) instead
  # of the pooled dispersion already shown on the Y axis. A tsv
  # regenerated after this column existed but without the aux stats
  # wired for it has the column present but blank ("" -> NA on read);
  # a tsv never regenerated at all (anything but entry-vs-anthology-
  # shinkokin as of 2026-10-15) doesn't have the column AT ALL -- d$
  # dispersion_difference is NULL, not a same-length NA vector, and
  # ifelse(is.na(NULL), 0, NULL) silently returns a ZERO-length vector,
  # not nrow(d) zeroes (confirmed 2026-10-15: broke entry-vs-gender-
  # first4.tsv with "replacement has 0 rows, data has 1655"). Both
  # cases -> uniform shade (0 everywhere), not an error.
  dispersion <- if (is.null(d$dispersion_difference)) {
    rep(0, nrow(d))
  } else {
    ifelse(is.na(d$dispersion_difference), 0, d$dispersion_difference)
  }
  # signed_range is both the real min/max for the legend's two end
  # labels AND now (2026-10-17, reverted off quantile ranking back to
  # plain linear normalization -- see .linear_ramp's own comment) the
  # range the color mapping itself is scaled against. abs_range stays
  # magnitude-only, kept for the per-type highlight ramps below
  # (.abs_ramp), which still encode "how different" via a single hue's
  # light/dark depth, not sign.
  signed_range <- range(dispersion)
  abs_range <- range(abs(dispersion))
  has_categories <- (!is.null(highlight) && nrow(highlight) > 0) ||
    (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) ||
    isTRUE(highlight_frontier)
  default_palette <- if (has_categories) default_palette else neutral_palette
  base_alpha <- if (has_categories) unselected_alpha else 0.85
  point_col <- adjustcolor(.linear_ramp(dispersion, signed_range, gradient_palette["lo"],
                                         gradient_palette["hi"]),
                            alpha.f = base_alpha)
  legend_n_steps <- 60
  legend_step_vals <- seq(signed_range[1], signed_range[2], length.out = legend_n_steps)
  legend_cols <- .linear_ramp(legend_step_vals, signed_range, gradient_palette["lo"],
                               gradient_palette["hi"])
  type_ramps <- list()
  if (!is.null(highlight) && nrow(highlight) > 0) {
    point_stem <- sub("【.*$", "", d$point)
    types <- unique(highlight$type)
    for (i in seq_along(types)) {
      type_name <- types[i]
      pal <- category_palettes[[((i - 1) %% length(category_palettes)) + 1]]
      terms <- highlight$term[highlight$type == type_name]
      is_match <- point_stem %in% terms
      if (!is.null(annotate_point_exclude)) {
        is_match <- is_match & !(d$point %in% annotate_point_exclude)
      }
      if (any(is_match)) {
        point_col[is_match] <- .abs_ramp(dispersion[is_match], abs_range, pal["light"], pal["dark"])
      }
      type_ramps[[type_name]] <- pal
    }
  }
  # 凡例の名前はCLI引数じゃなくtype列(テーブル側)から取る(2026-09-20確認)
  # ——highlight(単語版)のtype列と同じ扱い。type列が無い古い名単
  # (2列だけのTSV)にも後方互換で対応：その場合は単一カテゴリ"hl"として
  # 扱う。パレット色は単語版highlightで使った分の続きから割り当てる
  # (両方同時に使っても色が衝突しない)。
  if (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) {
    point_keys <- mapply(function(a, b) {
      paste(sort(sub("【.*$", "", c(a, b))), collapse = "\t")
    }, d$w1, d$w2)
    pair_types <- if (!is.null(highlight_pairs$type)) unique(highlight_pairs$type) else "hl"
    base_idx <- length(type_ramps)
    for (i in seq_along(pair_types)) {
      type_name <- pair_types[i]
      pal <- category_palettes[[((base_idx + i - 1) %% length(category_palettes)) + 1]]
      rows_i <- if (!is.null(highlight_pairs$type)) {
        highlight_pairs[highlight_pairs$type == type_name, ]
      } else {
        highlight_pairs
      }
      pair_keys <- unique(apply(rows_i[, c("w1", "w2")], 1,
                                 function(r) paste(sort(trimws(r)), collapse = "\t")))
      is_match <- point_keys %in% pair_keys
      if (any(is_match)) {
        point_col[is_match] <- .abs_ramp(dispersion[is_match], abs_range, pal["light"], pal["dark"])
      }
      type_ramps[[type_name]] <- pal
    }
  }
  frontier_match_idx <- NULL
  frontier_pal <- NULL
  if (isTRUE(highlight_frontier)) {
    frontier_match_idx <- .frontier_match_idx(d, frontier_n, frontier_candidate_idx,
                                               annotate_point_exclude)
    base_idx <- length(type_ramps)
    frontier_pal <- category_palettes[[(base_idx %% length(category_palettes)) + 1]]
    if (any(frontier_match_idx)) {
      point_col[frontier_match_idx] <- .abs_ramp(dispersion[frontier_match_idx], abs_range,
                                                  frontier_pal["light"], frontier_pal["dark"])
    }
    type_ramps[["frontier"]] <- frontier_pal
  }
  list(point_col = point_col, signed_range = signed_range, legend_cols = legend_cols,
       type_ramps = type_ramps, default_palette = default_palette,
       dispersion = dispersion, abs_range = abs_range,
       frontier_match_idx = frontier_match_idx, frontier_pal = frontier_pal,
       unselected_alpha = unselected_alpha)
}

# frontier色付け/frontierバッジ選定どちらにも使う"POOL_IDX（NULLなら全D）
# のうち、d$frequency降順でFRONTIER_N件"という判定を1箇所に集約——
# .compute_point_colors(着色用)と.resolve_inset_bundles(2026-11-07以降、
# インセットごとに同じ判定を自分のcand_idx上でやり直すため)の両方から
# 同じロジックを呼ぶ。返り値はnrow(d)長のlogical(どのdでも同じ行数を
# 前提にできるのは呼び出し側が常にPOOL_IDXを同じdへの添字として渡すため)。
.frontier_match_idx <- function(d, frontier_n, pool_idx, annotate_point_exclude) {
  pool_idx <- if (is.null(pool_idx)) seq_len(nrow(d)) else pool_idx
  ord <- pool_idx[order(-d$frequency[pool_idx])]
  top_n <- min(frontier_n, length(ord))
  top_idx <- ord[seq_len(top_n)]
  is_match <- seq_len(nrow(d)) %in% top_idx
  if (!is.null(annotate_point_exclude)) {
    is_match <- is_match & !(d$point %in% annotate_point_exclude)
  }
  is_match
}


.dedup_pair_rows <- function(rows) {
  if (is.null(rows$w1)) return(rows)
  has_assoc <- !is.na(rows$assoc_w1_to_w2) & !is.na(rows$assoc_w2_to_w1)
  rows <- rows[has_assoc, ]
  key <- paste(pmin(rows$w1, rows$w2), pmax(rows$w1, rows$w2), sep = "\t")
  rows[!duplicated(key), ]
}

.point_identity_key <- function(rows) {
  if (!is.null(rows$w1)) paste(pmin(rows$w1, rows$w2), pmax(rows$w1, rows$w2), sep = "\t")
  else rows$point
}

filter_to_annotated_points <- function(rows, badge_rows) {
  rows[.point_identity_key(rows) %in% .point_identity_key(badge_rows), ]
}

.pair_stem_key <- function(w1, w2) mapply(function(a, b) {
  paste(sort(sub("【.*$", "", c(a, b))), collapse = "\t")
}, w1, w2)

.select_badge_rows <- function(d, annotate_words, label_top_n, neutral,
                                annotate_nontarget_only, annotate_point_exclude,
                                badge_selection = "words", frontier_n = 15,
                                annotate_pairs = NULL) {
  if (identical(badge_selection, "frontier")) {
    matched <- .dedup_pair_rows(d)
    matched <- matched[order(-matched$frequency), ]
    if (nrow(matched) > frontier_n) {
      matched <- matched[seq_len(frontier_n), ]
    }
    return(matched)
  }
  if (identical(badge_selection, "pairs")) {
    stopifnot(!is.null(d$w1), !is.null(annotate_pairs))
    target_keys <- .pair_stem_key(annotate_pairs$w1, annotate_pairs$w2)
    matched <- d[.pair_stem_key(d$w1, d$w2) %in% target_keys, ]
    return(.dedup_pair_rows(matched))
  }
  if (identical(badge_selection, "points")) {
    stopifnot(!is.null(annotate_words))
    matched <- d[d$point %in% annotate_words, ]
    hi <- matched[matched$keyness > neutral, ]
    lo <- matched[matched$keyness <= neutral, ]
    return(.dedup_pair_rows(rbind(hi, lo)))
  }
  if (!is.null(annotate_words)) {
    point_stem <- sub("【.*$", "", d$point)
    matched <- d[point_stem %in% annotate_words, ]
    if (annotate_nontarget_only) matched <- matched[matched$n_target == 0, ]
    if (!is.null(annotate_point_exclude)) {
      matched <- matched[!(matched$point %in% annotate_point_exclude), ]
    }
    hi <- matched[matched$keyness > neutral, ]
    lo <- matched[matched$keyness <= neutral, ]
  } else {
    # Ranks by d$frequency, not d$keyness -- since the axis swap in
    # draw_point_vs_target() (always applied before this function ever
    # runs; 2026-11-05: d$frequency <- d$keyness, d$keyness <- the
    # adjusted frequency (the original node_size column), d$node_size <-
    # d$dispersion), d$frequency is the one still holding the original
    # forward KLD (point->target), while d$keyness now holds adjusted
    # frequency (non-negative, not a KLD at all) instead of dispersion.
    # Default labels pick the forward-KLD extremes (most positive/negative
    # point->target divergence), not frequency/dispersion extremes,
    # confirmed with the user 2026-10-15 (axis choice) and 2026-11-05
    # (frequency/dispersion swap, see the file header comment).
    hi <- d[order(-d$frequency), ][seq_len(min(label_top_n, nrow(d))), ]
    lo <- d[order(d$frequency), ][seq_len(min(label_top_n, nrow(d))), ]
  }
  .dedup_pair_rows(rbind(hi, lo))
}

.aggregate_locator_groups <- function(rows) {
  if (nrow(rows) == 0) return(rows)
  is_pair <- !is.null(rows$w1)
  kind_key <- if (is_pair) {
    if (!is.null(rows$field1)) paste(rows$field1, rows$field2, sep = "|") else ""
  } else {
    if (!is.null(rows$field)) rows$field else ""
  }
  key <- paste(rows$n_target, rows$n_not_target, kind_key, sep = "\t")
  idx_by_key <- split(seq_len(nrow(rows)), key)
  pair_member_cols <- c("w1", "w2", "freq1", "freq2", "assoc_w1_to_w2", "assoc_w2_to_w1")
  if (is_pair && !is.null(rows$field1)) pair_member_cols <- c(pair_member_cols, "field1", "field2")
  out_rows <- lapply(idx_by_key, function(idx) {
    rep_row <- rows[idx[1], , drop = FALSE]
    rep_row$group_size <- length(idx)
    rep_row$members <- if (is_pair) {
      list(rows[idx, pair_member_cols])
    } else {
      list(rows$point[idx])
    }
    rep_row
  })
  do.call(rbind, out_rows)
}

.badge_color_from_point <- function(rows) {
  if (nrow(rows) == 0) {
    return(list(fill = character(0), text_col = character(0), border_col = character(0),
                 word_col = character(0)))
  }
  fill <- rows$point_col
  rgba <- col2rgb(fill, alpha = TRUE)
  a <- rgba["alpha", ] / 255
  r_eff <- rgba["red", ] * a + 255 * (1 - a)
  g_eff <- rgba["green", ] * a + 255 * (1 - a)
  b_eff <- rgba["blue", ] * a + 255 * (1 - a)
  luminance <- (0.299 * r_eff + 0.587 * g_eff + 0.114 * b_eff) / 255
  dark <- luminance < 0.45
  word_col <- rgb(rgba["red", ], rgba["green", ], rgba["blue", ], maxColorValue = 255)
  list(fill = fill, word_col = word_col,
       text_col = ifelse(dark, "white", "#111111"),
       border_col = ifelse(dark, "white", "#555555"))
}


.render_pair_network_svg <- function(members, fontsize = 11, network_scale = 0.01, network_node_gap = 1,
                                       filter_edge_count = 500, top_pct = 0.01) {
  if (nrow(members) > filter_edge_count) {
    strength <- pmax(abs(members$assoc_w1_to_w2), abs(members$assoc_w2_to_w1))
    thresh <- quantile(strength, 1 - top_pct)
    members <- members[strength >= thresh, ]
  }
  field1 <- if (!is.null(members$field1)) members$field1 else NA_character_
  field2 <- if (!is.null(members$field2)) members$field2 else NA_character_
  members$w1 <- .concept_display_text(members$w1, field1)
  members$w2 <- .concept_display_text(members$w2, field2)
  nodes <- unique(c(members$w1, members$w2))
  esc <- function(s) gsub('"', '\\\\"', s, fixed = TRUE)
  max_abs <- max(abs(members$assoc_w1_to_w2), abs(members$assoc_w2_to_w1))
  lw <- function(v) sprintf("%.2f", pmax(0.8, pmin(4, abs(v) / max_abs * 4)))
  node_col <- setNames(rep("#111111", length(nodes)), nodes)
  edge_decls <- character(0)
  for (i in seq_len(nrow(members))) {
    w1 <- members$w1[i]; w2 <- members$w2[i]
    v12 <- members$assoc_w1_to_w2[i]; v21 <- members$assoc_w2_to_w1[i]
    col1 <- if (v12 > v21) "#111111" else "#CC2200"
    col2 <- if (v12 > v21) "#CC2200" else "#111111"
    if (col1 == "#CC2200") node_col[[w1]] <- "#CC2200"
    if (col2 == "#CC2200") node_col[[w2]] <- "#CC2200"
    edge_decls <- c(edge_decls,
      sprintf('  "%s" -> "%s" [color="%s", penwidth=%s];', esc(w1), esc(w2), col1, lw(v12)),
      sprintf('  "%s" -> "%s" [color="%s", penwidth=%s];', esc(w2), esc(w1), col2, lw(v21)))
  }
  node_field <- setNames(rep(NA_character_, length(nodes)), nodes)
  if (!is.null(members$field1)) {
    for (i in seq_len(nrow(members))) {
      if (is.na(node_field[[members$w1[i]]])) node_field[[members$w1[i]]] <- members$field1[i]
      if (is.na(node_field[[members$w2[i]]])) node_field[[members$w2[i]]] <- members$field2[i]
    }
  }
  node_shape <- vapply(node_field[nodes], function(f) {
    if (is.na(f)) return("ellipse")
    switch(f, entry = "ellipse", concept = "box", topic = "diamond", "ellipse")
  }, character(1))
  node_decls <- sprintf('  "%s" [color="%s", shape=%s];', esc(nodes), node_col[nodes], node_shape)
  margin_w <- 0.06 * network_scale; margin_h <- 0.03 * network_scale
  dot_src <- c("digraph g {",
               "  layout=neato;", "  overlap=scale;",
               sprintf('  sep="+%g";', network_node_gap), "  bgcolor=transparent;",
               sprintf('  node [shape=ellipse, style=filled, fillcolor=white, fontname="Hiragino Mincho ProN", fontsize=%d, penwidth=1.6, width=0.15, height=0.1, margin="%g,%g"];',
                       fontsize, margin_w, margin_h),
               '  edge [fontname="Hiragino Mincho ProN"];',
               "  splines=curved;",
               node_decls, edge_decls, "}")
  svg_lines <- system2("neato", args = "-Tsvg", input = dot_src, stdout = TRUE)
  status <- attr(svg_lines, "status")
  if (!is.null(status) && status != 0) {
    stop("neato -Tsvg failed (需要pkgs.graphviz在PATH里，见flake.nix)")
  }
  paste(svg_lines, collapse = "\n")
}

.parse_svg_dims_pt <- function(svg_txt) {
  header <- regmatches(svg_txt, regexpr('width="[0-9.]+pt" height="[0-9.]+pt"', svg_txt))
  nums <- as.numeric(regmatches(header, gregexpr("[0-9.]+", header))[[1]])
  list(w = nums[1], h = nums[2])
}

.concept_display_labels_path <- function() {
  Sys.getenv("CONCEPT_DISPLAY_LABELS", "data/processed/concept-display-labels.tsv")
}

.concept_display_label_map <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      path <- .concept_display_labels_path()
      cache <<- if (file.exists(path)) {
        tbl <- read.table(path, sep = "\t", header = FALSE, quote = "",
                           stringsAsFactors = FALSE, encoding = "UTF-8",
                           col.names = c("key", "display"))
        setNames(tbl$display, tbl$key)
      } else {
        character(0)
      }
    }
    cache
  }
})

# TEI本体でXML NameCharの範囲外の文字(例えば"*"=U+002A、"+"=U+002B)は
# xml:id内で"_x<4桁hex>_"形式にエスケープされたまま、measured_value/w1/
# w2までそのまま素通りしてくる(entry2l2entry_exclusion.pyの既存コメント・
# 実データ"の_x002A_【の】"参照) -- 表示の最後の一歩でだけ元の文字に戻す
# (pair keyを"+"で分割するコードはこの関数を通す前の、まだエスケープされた
# ままの文字列を読む -- 通した後の文字列を再度"+"で分割してはいけない、
# "+"自体がエスケープを解いた結果として出てくることがあるため)。
.unescape_xml_id_chars <- function(s) {
  vapply(s, function(x) {
    if (is.na(x) || !grepl("_x[0-9A-Fa-f]{4}_", x)) return(x)
    matches <- unique(regmatches(x, gregexpr("_x[0-9A-Fa-f]{4}_", x))[[1]])
    for (mt in matches) {
      ch <- intToUtf8(strtoi(substr(mt, 3, 6), base = 16L))
      x <- gsub(mt, ch, x, fixed = TRUE)
    }
    x
  }, character(1), USE.NAMES = FALSE)
}

# 部立大分類(budatemajor、field_label_by_type()の"topic")の表示は【】では
# なく『』で囲む -- 元データのbudatemajorラベル自体には括弧が一切ない
# (例"その他")、entry側の既存の【漢字】記法と視覚的に区別するための表示
#専用の追加(2026-10-18ユーザー要望)。FIELD_TYPEがこの文字列以外(NA含む)
# の時は何もしない -- 複合済みの"w1・w2"文字列をもう一度この関数に通す
# 呼び出し(後述の"word"フォールバック)はFIELD_TYPEにNAを渡し、二重に
# 括弧で囲まれることを防ぐ。
.wrap_topic_label <- function(text, field_type) {
  is_topic <- !is.na(field_type) & field_type == "部立大分類"
  ifelse(is_topic & !is.na(text), paste0("『", text, "』"), text)
}

.pair_display_sep <- function(field1, field2) {
  is_budate <- (!is.na(field1) & field1 == "部立大分類") | (!is.na(field2) & field2 == "部立大分類")
  ifelse(is_budate, "，", "+")
}

.concept_display_text <- function(s, field_type = NA_character_) {
  hit <- .concept_display_label_map()[s]
  text <- ifelse(!is.na(hit), hit, sub("（BG-[^）]*）", "", s))
  text <- .unescape_xml_id_chars(text)
  text <- .wrap_topic_label(text, field_type)
  text
}

.format_display_columns <- function(d) {
  if (!is.null(d$w1)) {
    field1 <- if (!is.null(d$field1)) d$field1 else NA_character_
    field2 <- if (!is.null(d$field2)) d$field2 else NA_character_
    d$w1 <- .concept_display_text(d$w1, field1)
    d$w2 <- .concept_display_text(d$w2, field2)
  } else if (!is.null(d$point)) {
    d$point <- .concept_display_text(d$point, d$field)
  }
  d
}

.rbind_fill <- function(dfs) {
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) df[[col]] <- NA
    df[all_cols]
  })
  do.call(rbind, dfs2)
}

.badge_geometry <- function(rows, value_label, not_value_label, badge_cex, fam, label_gap_in,
                              max_group_members = 50) {
  is_pair <- !is.null(rows$w1)
  has_groups <- !is.null(rows$group_size)
  lapply(seq_len(nrow(rows)), function(i) {
    count_label <- sprintf("%s:%d / %s:%d", value_label, rows$n_target[i],
                            not_value_label, rows$n_not_target[i])
    count_w_in <- strwidth(count_label, cex = badge_cex, family = fam, units = "inches")
    count_h_in <- strheight(count_label, cex = badge_cex, family = fam, units = "inches")
    box_w_in <- count_w_in + count_w_in * 0.35
    box_h_in <- count_h_in + count_h_in * 1.2
    group_size <- if (has_groups) rows$group_size[i] else 1

    if (is_pair) {
      pair_cols <- intersect(c("w1", "w2", "field1", "field2"), names(rows))
      members_i <- if (has_groups) {
        rows$members[[i]]
      } else {
        rows[i, pair_cols, drop = FALSE]
      }
      field1 <- if (!is.null(members_i$field1)) members_i$field1 else NA_character_
      field2 <- if (!is.null(members_i$field2)) members_i$field2 else NA_character_
      members <- paste0(.concept_display_text(members_i$w1, field1),
                         .pair_display_sep(field1, field2),
                         .concept_display_text(members_i$w2, field2))
    } else if (group_size > 1) {
      members <- .concept_display_text(rows$members[[i]], rows$field[i])
    } else {
      members <- NULL
    }

    if (!is.null(members) && length(members) > 1) {
      shown <- if (length(members) > max_group_members) {
        c(members[seq_len(max_group_members)],
          sprintf("…計%d件", length(members)))
      } else {
        members
      }
      line_w_in <- vapply(shown, strwidth, numeric(1), cex = badge_cex, family = fam, units = "inches")
      line_h_in <- vapply(shown, strheight, numeric(1), cex = badge_cex, family = fam, units = "inches")
      line_gap_in <- label_gap_in * 1.8
      return(list(kind = "word_list", word_list = shown, count_label = count_label,
                   box_w_in = box_w_in, box_h_in = box_h_in,
                   combined_w_in = max(max(line_w_in), box_w_in),
                   combined_h_in = box_h_in + label_gap_in + sum(line_h_in) + line_gap_in * length(shown)))
    }
    # ここに来る時点でis_pairなら(group_size==1でも)membersは常に非NULLで
    # 既に.concept_display_text済み("w1・w2"結合後の文字列) -- 二重に
    # 『』を付けないようFIELD_TYPEはNAを渡す(.wrap_topic_labelはNAなら
    # 素通り)。!is_pairならここに来る時点でgroup_size==1が保証されている
    # (group_size>1は上のword_list分岐がreturn済み)ので、membersはNULL、
    # rows$point[i]を生のまま.field[i]付きで処理する。
    word_field <- if (is_pair) NA_character_ else rows$field[i]
    word <- .concept_display_text(if (!is.null(members)) members[1] else rows$point[i], word_field)
    word_w_in <- strwidth(word, cex = badge_cex, family = fam, units = "inches")
    word_h_in <- strheight(word, cex = badge_cex, family = fam, units = "inches")
    list(kind = "word", word = word, count_label = count_label,
         box_w_in = box_w_in, box_h_in = box_h_in,
         combined_w_in = max(word_w_in, box_w_in),
         combined_h_in = box_h_in + label_gap_in + word_h_in)
  })
}


.layout_badges_repel <- function(combined_w, combined_h, px_in, py_in, pin, badge_gap_in) {
  n <- length(combined_w)
  x1 <- px_in - 0.01 * pin[1] - combined_w
  y1 <- py_in + badge_gap_in
  box_pad_y_in <- badge_gap_in * 1.5
  boxes <- cbind(x1 = x1, y1 = y1 - box_pad_y_in,
                 x2 = x1 + combined_w, y2 = y1 + combined_h + box_pad_y_in)
  repelled <- ggrepel:::repel_boxes2(
    data_points = cbind(px_in, py_in), point_size = rep(0.02, n),
    point_padding_x = badge_gap_in, point_padding_y = badge_gap_in,
    boxes = boxes, xlim = c(0, pin[1]), ylim = c(0, pin[2]),
    hjust = rep(0.5, n), vjust = rep(0.5, n),
    force_push = 1e-4, force_pull = 1e-6, max_time = 3,
    max_overlaps = 20, max_iter = 3000, direction = "both"
  )
  list(x1 = repelled$x, y1 = repelled$y + box_pad_y_in)
}

.largest_empty_rect <- function(occupied) {
  nr <- nrow(occupied); nc <- ncol(occupied)
  heights <- integer(nc)
  best <- list(area = 0, r1 = 1, r2 = 0, c1 = 1, c2 = 0)
  for (r in seq_len(nr)) {
    heights <- ifelse(occupied[r, ], 0L, heights + 1L)
    hh <- c(heights, 0L)
    stack <- integer(0)
    for (c in seq_len(nc + 1)) {
      while (length(stack) > 0 && hh[stack[length(stack)]] >= hh[c]) {
        top <- stack[length(stack)]
        stack <- stack[-length(stack)]
        width <- if (length(stack) == 0) c - 1 else c - stack[length(stack)] - 1
        area <- hh[top] * width
        if (area > best$area) {
          best <- list(area = area, r1 = r - hh[top] + 1, r2 = r,
                       c1 = if (length(stack) == 0) 1 else stack[length(stack)] + 1,
                       c2 = c - 1)
        }
      }
      stack <- c(stack, c)
    }
  }
  best
}

.layout_badges_empty_region <- function(combined_w, combined_h, row_order,
                                          px_in, py_in,
                                          all_px_in, all_py_in, point_radius_in,
                                          pin, badge_gap_in, res = 60) {
  n <- length(combined_w)
  cell_w <- pin[1] / res
  cell_h <- pin[2] / res
  occupied <- matrix(FALSE, nrow = res, ncol = res)
  cx_idx <- pmin(res, pmax(1, ceiling(all_px_in / cell_w)))
  cy_idx <- pmin(res, pmax(1, ceiling((pin[2] - all_py_in) / cell_h)))
  pad_x <- max(1, ceiling(point_radius_in / cell_w))
  pad_y <- max(1, ceiling(point_radius_in / cell_h))
  for (k in seq_along(cx_idx)) {
    xr <- max(1, cx_idx[k] - pad_x):min(res, cx_idx[k] + pad_x)
    yr <- max(1, cy_idx[k] - pad_y):min(res, cy_idx[k] + pad_y)
    occupied[yr, xr] <- TRUE
  }
  # text_scale>1(2026-11-14)：badgeが大きくなり空白が減るため、凡例の長い
  # 種別名(例「部立×語→女性歌KLD上位10」、約1.4in)の上にまでbadgeが
  # 置かれるようになった——予約幅を凡例の実寸まで広げる。既定1では従来どおり。
  legend_w_in <- if (.text_scale() > 1) 1.45 else 0.75
  legend_h_in <- if (.text_scale() > 1) 2.5 else 2.2
  lx1 <- max(1, ceiling((pin[1] - legend_w_in) / cell_w))
  ly2 <- max(1, ceiling(legend_h_in / cell_h))
  occupied[1:min(res, ly2), lx1:res] <- TRUE

  x1v <- numeric(n); y1v <- numeric(n)
  for (i in row_order) {
    best <- .largest_empty_rect(occupied)
    region <- c((best$c1 - 1) * cell_w, pin[2] - best$r2 * cell_h,
                best$c2 * cell_w, pin[2] - (best$r1 - 1) * cell_h)
    w <- combined_w[i] + badge_gap_in
    h <- combined_h[i] + badge_gap_in
    x1 <- max(region[1], min(px_in[i] - w / 2, region[3] - w))
    y1 <- max(region[2], min(py_in[i] - h / 2, region[4] - h))
    x1 <- max(0, min(x1, pin[1] - w))
    y1 <- max(0, min(y1, pin[2] - h))
    x1v[i] <- x1
    y1v[i] <- y1

    occ_c1 <- max(1, floor(x1 / cell_w) + 1)
    occ_c2 <- min(res, ceiling((x1 + w) / cell_w))
    occ_r1 <- max(1, floor((pin[2] - (y1 + h)) / cell_h) + 1)
    occ_r2 <- min(res, ceiling((pin[2] - y1) / cell_h))
    occupied[occ_r1:occ_r2, occ_c1:occ_c2] <- TRUE
  }
  list(x1 = x1v, y1 = y1v)
}


.draw_connector_line <- function(px_in, py_in, combined_x1_in, combined_y1_in,
                                  combined_w, combined_h, box_w_in, box_h_in, col,
                                  in_to_data_x, in_to_data_y) {
  combined_x2_in <- combined_x1_in + combined_w
  combined_y2_in <- combined_y1_in + combined_h
  box_cx_in <- (combined_x1_in + combined_x2_in) / 2
  box_cy_in <- (combined_y1_in + combined_y2_in) / 2
  box_x1_in <- box_cx_in - box_w_in / 2
  box_x2_in <- box_cx_in + box_w_in / 2
  enters_vertically <- abs(py_in - box_cy_in) >= abs(px_in - box_cx_in)
  if (enters_vertically && py_in < box_cy_in) enters_vertically <- FALSE
  if (enters_vertically) {
    anchor_x_in <- box_cx_in
    anchor_y_in <- combined_y2_in
  } else {
    anchor_y_in <- combined_y2_in - 0.25 * box_h_in
    anchor_x_in <- if (px_in <= box_cx_in) box_x1_in else box_x2_in
  }
  dx_in <- anchor_x_in - px_in
  dy_in <- anchor_y_in - py_in
  bend_frac <- 0.6
  if (enters_vertically) {
    mid_x_in <- anchor_x_in
    mid_y_in <- py_in + dy_in * bend_frac
  } else {
    mid_y_in <- anchor_y_in
    mid_x_in <- px_in + dx_in * bend_frac
  }
  lines(in_to_data_x(c(px_in, mid_x_in, anchor_x_in)),
        in_to_data_y(c(py_in, mid_y_in, anchor_y_in)),
        col = col, lwd = 0.8, xpd = NA)
}

.draw_badge_box <- function(g, combined_x1_in, combined_y1_in, combined_w, combined_h,
                             fill, border, text_col, word_col, badge_cex, fam,
                             label_gap_in, x_per_in, y_per_in, in_to_data_x, in_to_data_y) {
  cx_in <- combined_x1_in + combined_w / 2
  cx <- in_to_data_x(cx_in)
  count_top_in <- combined_y1_in + combined_h
  count_bot_in <- count_top_in - g$box_h_in
  y_top <- in_to_data_y(count_top_in)
  box_w <- g$box_w_in * x_per_in
  box_h <- g$box_h_in * y_per_in

  rect(cx - box_w / 2, y_top - box_h, cx + box_w / 2, y_top,
       col = fill, border = border, xpd = NA)
  text(cx, y_top - box_h / 2, labels = g$count_label, cex = badge_cex,
       family = fam, col = text_col, xpd = NA)
  # Word text itself stays fixed black, not word_col (the point's own
  # opaque color -- still used for the connector line below, unchanged)
  # -- badge selection no longer tracks the color ramp's dark end (it's
  # picked by forward-KLD extremes now, see .select_badge_rows), so a
  # labeled word can easily land near abs(dispersion_difference)==0,
  # the ramp's lightest/near-white end, which is unreadable as text.
  # Confirmed with the user 2026-10-15.
  word_text_col <- "#111111"
  if (identical(g$kind, "word_list")) {
    line_gap_in <- label_gap_in * 1.8
    y_in <- count_bot_in - label_gap_in
    for (line in g$word_list) {
      h_in <- strheight(line, cex = badge_cex, family = fam, units = "inches")
      text(cx, in_to_data_y(y_in), labels = line,
           cex = badge_cex, family = fam, col = word_text_col, adj = c(0.5, 1), xpd = NA)
      y_in <- y_in - h_in - line_gap_in
    }
  } else {
    text(cx, in_to_data_y(count_bot_in - label_gap_in), labels = g$word,
         cex = badge_cex, family = fam, col = word_text_col, adj = c(0.5, 1), xpd = NA)
  }
}

.draw_badges <- function(rows, fill, border, text_col, word_col, badge_layout,
                          d_all, value_label, not_value_label, badge_cex, fam,
                          label_gap_in, x_per_in, y_per_in, data_to_in_x, data_to_in_y,
                          in_to_data_x, in_to_data_y, pin, neutral, max_group_members = 50,
                          tie_boxes = list()) {
  n <- nrow(rows)
  if (n == 0) return(invisible(NULL))
  if (length(fill) == 1) fill <- rep(fill, n)
  if (length(border) == 1) border <- rep(border, n)
  if (length(text_col) == 1) text_col <- rep(text_col, n)
  if (length(word_col) == 1) word_col <- rep(word_col, n)

  geoms <- .badge_geometry(rows, value_label, not_value_label, badge_cex, fam, label_gap_in,
                            max_group_members = max_group_members)
  combined_w <- vapply(geoms, `[[`, numeric(1), "combined_w_in")
  combined_h <- vapply(geoms, `[[`, numeric(1), "combined_h_in")
  box_w <- vapply(geoms, `[[`, numeric(1), "box_w_in")
  box_h <- vapply(geoms, `[[`, numeric(1), "box_h_in")
  px_in <- data_to_in_x(rows$frequency)
  py_in <- data_to_in_y(rows$keyness)
  badge_gap_in <- label_gap_in
  row_order <- order(-abs(rows$keyness - neutral))

  resolved <- if (identical(badge_layout, "repel")) {
    .layout_badges_repel(combined_w, combined_h, px_in, py_in, pin, badge_gap_in)
  } else {
    all_px_in <- data_to_in_x(d_all$frequency)
    all_py_in <- data_to_in_y(d_all$keyness)
    point_radius_in <- 0.04
    .layout_badges_empty_region(combined_w, combined_h, row_order, px_in, py_in,
                                 all_px_in, all_py_in, point_radius_in, pin, badge_gap_in)
  }

  for (i in seq_len(n)) {
    .draw_badge_box(geoms[[i]], resolved$x1[i], resolved$y1[i], combined_w[i], combined_h[i],
                     fill[i], border[i], text_col[i], word_col[i], badge_cex, fam,
                     label_gap_in, x_per_in, y_per_in, in_to_data_x, in_to_data_y)
  }
  # ROWSの各点がtie_boxes(.compute_x_tie_boxesの出力、{key: {x1,y1,x2,
  # y2}})のどれかに属していれば、連結線の出発点はその点自身ではなく箱の
  # 辺にする(2026-10-17、ユーザー確認済み：「線がもし框なら、框のedgeから
  # 出発。もし点なら、そのまま」)。出発点を「箱の中でbadge自身の位置に
  # 一番近い点」にクランプするだけで済む -- badgeは普通箱の外にあるので、
  # x/yそれぞれ独立に[box_x1,box_x2]/[box_y1,box_y2]へクランプした結果は
  # そのまま箱の辺(角を含む)に乗る、標準的な「矩形上の最近点」の求め方。
  start_x_in <- px_in
  start_y_in <- py_in
  if (length(tie_boxes) > 0) {
    keys <- paste(rows$n_target, rows$n_not_target, sep = "/")
    for (i in seq_len(n)) {
      box <- tie_boxes[[keys[i]]]
      if (is.null(box)) next
      box_x_in <- data_to_in_x(c(box$x1, box$x2))
      box_y_in <- data_to_in_y(c(box$y1, box$y2))
      badge_cx_in <- resolved$x1[i] + combined_w[i] / 2
      badge_cy_in <- resolved$y1[i] + combined_h[i] / 2
      start_x_in[i] <- min(max(badge_cx_in, min(box_x_in)), max(box_x_in))
      start_y_in[i] <- min(max(badge_cy_in, min(box_y_in)), max(box_y_in))
    }
  }
  for (i in seq_len(n)) {
    .draw_connector_line(start_x_in[i], start_y_in[i], resolved$x1[i], resolved$y1[i],
                          combined_w[i], combined_h[i], box_w[i], box_h[i], word_col[i],
                          in_to_data_x, in_to_data_y)
  }
  invisible(NULL)
}


.draw_svg_calibration_marks <- function(in_to_data_x, in_to_data_y, pin) {
  mark_col <- "#010203"
  dx_per_in <- in_to_data_x(1) - in_to_data_x(0)
  dy_per_in <- in_to_data_y(1) - in_to_data_y(0)
  mark_half_in <- 0.001
  marks_in <- list(c(0.15, 0.15), c(max(pin[1] - 0.15, 0.15), max(pin[2] - 0.15, 0.15)))
  for (p in marks_in) {
    x0 <- in_to_data_x(p[1]); y0 <- in_to_data_y(p[2])
    dx <- mark_half_in * dx_per_in; dy <- mark_half_in * dy_per_in
    rect(x0 - dx, y0 - dy, x0 + dx, y0 + dy, col = mark_col, border = NA, xpd = NA)
  }
  marks_in
}

.splice_network_glyphs <- function(output, pending, marks_in) {
  if (length(pending) == 0) return(invisible(NULL))
  if (!grepl("\\.svg$", output, ignore.case = TRUE)) {
    stop("badge_layout=\"network\" needs svg output (graphviz glyphs are ",
         "spliced in as SVG XML text) -- png output doesn't support this badge layout")
  }
  svg_txt <- paste(readLines(output, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  lines_all <- strsplit(svg_txt, "\n")[[1]]
  mark_lines <- grep("0.392157%, 0.784314%, 1.176471%", lines_all, fixed = TRUE, value = TRUE)
  expected_marks <- 2 * length(pending)
  if (length(mark_lines) != expected_marks) {
    stop(sprintf("找不到svg校准标记(应该有%d个，每个网络面板自己2个)，网络badge拼接失败——可能是svg输出格式变了",
                 expected_marks))
  }
  get_center <- function(line) {
    d_str <- sub('.*d="([^"]*)".*', "\\1", line)
    nums <- as.numeric(regmatches(d_str, gregexpr("[0-9]+\\.?[0-9]*", d_str))[[1]])
    c(mean(nums[c(1, 3)]), mean(nums[c(2, 6)]))
  }
  lines_all <- lines_all[!grepl("0.392157%, 0.784314%, 1.176471%", lines_all, fixed = TRUE)]

  all_blocks <- character(0)
  for (i in seq_along(pending)) {
    raw1 <- get_center(mark_lines[2 * i - 1]); raw2 <- get_center(mark_lines[2 * i])
    in1 <- marks_in[[i]][[1]]; in2 <- marks_in[[i]][[2]]
    scale_x <- (raw2[1] - raw1[1]) / (in2[1] - in1[1])
    scale_y <- (raw2[2] - raw1[2]) / (in2[2] - in1[2])
    if (abs(scale_x - 72) > 0.5 || abs(scale_y + 72) > 0.5) {
      stop(sprintf("面板%d的svg校准标记反解出的scale(%.3f, %.3f)偏离预期的(72, -72)太多，网络badge拼接的坐标换算前提不成立，需要重新核对",
                   i, scale_x, scale_y))
    }
    offset_x <- raw1[1] - scale_x * in1[1]
    offset_y <- raw1[2] - scale_y * in1[2]
    affine_x <- function(ix) scale_x * ix + offset_x
    affine_y <- function(iy) scale_y * iy + offset_y

    blocks <- vapply(pending[[i]], function(p) {
      x1_raw <- affine_x(p$x1_in)
      y_top_raw <- affine_y(p$y1_in + p$h_in)
      w_raw <- p$w_in * 72
      h_raw <- p$h_in * 72
      inner <- sub("^.*?(<g )", "\\1", p$svg_txt)
      inner <- sub("</svg>\\s*$", "", inner)
      sprintf('<svg x="%f" y="%f" width="%f" height="%f" viewBox="0 0 %f %f">\n%s\n</svg>',
              x1_raw, y_top_raw, w_raw, h_raw, p$gv_w_pt, p$gv_h_pt, inner)
    }, character(1))
    all_blocks <- c(all_blocks, blocks)
  }

  svg_txt <- paste(lines_all, collapse = "\n")
  svg_txt <- sub("</svg>\\s*$", paste0(paste(all_blocks, collapse = "\n"), "\n</svg>\n"), svg_txt)
  writeLines(svg_txt, output, useBytes = TRUE)
  invisible(NULL)
}


.UI_LABELS <- list(
  ja = c(frequency = "総頻度", dispersion = "分散度", dispersion_range = "分散度差",
         frontier = "対象"),
  en = c(frequency = "total frequency", dispersion = "dispersion", dispersion_range = "dispersion range",
         frontier = "target")
)
.ui_label <- function(key, label_lang = "ja") {
  tbl <- .UI_LABELS[[label_lang]]
  if (is.null(tbl)) tbl <- .UI_LABELS[["ja"]]
  unname(tbl[key])
}

.draw_legend <- function(type_ramps, signed_range, legend_cols, default_label, pin,
                          in_to_data_x, in_to_data_y, default_palette = default_color_palette(),
                          label_lang = "ja") {
  legend_w_in <- 0.11
  legend_h_in <- 0.7
  legend_x2_in <- pin[1] - 0.05
  legend_x1_in <- legend_x2_in - legend_w_in
  legend_y1_in <- pin[2] - legend_h_in - 1.3
  n_steps <- length(legend_cols)
  step_h_in <- legend_h_in / n_steps
  for (i in seq_len(n_steps)) {
    rect(in_to_data_x(legend_x1_in), in_to_data_y(legend_y1_in + (i - 1) * step_h_in),
         in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + i * step_h_in),
         col = legend_cols[i], border = NA, xpd = NA)
  }
  rect(in_to_data_x(legend_x1_in), in_to_data_y(legend_y1_in),
       in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + legend_h_in),
       border = "#888888", xpd = NA)
  legend_cex <- 0.5
  text(in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + legend_h_in),
       .ui_label("dispersion_range", label_lang), cex = legend_cex, adj = c(1, -0.6), xpd = NA)
  text(in_to_data_x(legend_x1_in - 0.03), in_to_data_y(legend_y1_in + legend_h_in),
       format(signed_range[2], digits = 3), cex = legend_cex, adj = c(1, 1), xpd = NA)
  text(in_to_data_x(legend_x1_in - 0.03), in_to_data_y(legend_y1_in),
       format(signed_range[1], digits = 3), cex = legend_cex, adj = c(1, 0), xpd = NA)

  swatch_size_in <- 0.06
  swatch_gap_in <- 0.03
  if (length(type_ramps) == 0) return(invisible(NULL))
  all_labels <- c(list(default_palette), type_ramps)
  names(all_labels)[1] <- default_label
  label_right_in <- pin[1] - 0.05
  for (i in seq_along(all_labels)) {
    y_in <- legend_y1_in - swatch_gap_in - (i - 1) * (swatch_size_in + swatch_gap_in) - swatch_size_in
    pal <- all_labels[[i]]
    swatch_col <- colorRampPalette(c(pal["light"], pal["dark"]))(2)[2]
    label <- names(all_labels)[i]
    if (identical(label, "frontier")) label <- .ui_label("frontier", label_lang)
    text(in_to_data_x(label_right_in), in_to_data_y(y_in + swatch_size_in / 2),
         label, cex = legend_cex, adj = c(1, 0.5), xpd = NA)
    label_w_in <- strwidth(label, cex = legend_cex, units = "inches")
    swatch_x2_in <- label_right_in - label_w_in - swatch_gap_in
    swatch_x1_in <- swatch_x2_in - swatch_size_in
    rect(in_to_data_x(swatch_x1_in), in_to_data_y(y_in),
         in_to_data_x(swatch_x2_in), in_to_data_y(y_in + swatch_size_in),
         col = swatch_col, border = NA, xpd = NA)
  }
}

# 凡例は代表値3個(最小/中央値/最大)。中央値は偶数個のとき2つの中間値の
# 平均になり、実データのどの行の値とも一致しないことがある -- 「一番近い実際の行」を探してその行
# 自身の値を読む(2026-10-17確認済みの「表生成時に1回だけ計算、Rで再計算
# しない」という原則は変わらず、dispersion自体もtsv生成時に書き込まれた
# 列のまま)。
.draw_size_legend <- function(size_value, pin, in_to_data_x, in_to_data_y, pch = 1,
                               label_lang = "ja") {
  legend_cex <- 0.5
  x_in <- pin[1] - 0.15
  unique_vals <- sort(unique(size_value))
  target_values <- if (length(unique_vals) >= 3) {
    c(unique_vals[1], unique_vals[ceiling(length(unique_vals) / 2)],
      unique_vals[length(unique_vals)])
  } else {
    unique_vals
  }
  closest_idx <- vapply(target_values, function(v) which.min(abs(size_value - v)), integer(1))
  sizes <- size_value[closest_idx]
  labels <- sprintf("%.3f", sizes)
  gap_in <- 0.16
  y_top_in <- pin[2] - 0.22
  text(in_to_data_x(pin[1] - 0.05), in_to_data_y(pin[2] - 0.06),
       .ui_label("dispersion", label_lang), cex = legend_cex, adj = c(1, 0), xpd = NA)
  for (i in seq_along(labels)) {
    y_in <- y_top_in - (i - 1) * gap_in
    points(in_to_data_x(x_in), in_to_data_y(y_in), pch = pch,
           cex = .size_scale(sizes[i]), col = "#666666", xpd = NA)
    text(in_to_data_x(x_in - 0.12), in_to_data_y(y_in), labels[i],
         cex = legend_cex, adj = c(1, 0.5), xpd = NA)
  }
}

.representative_pch <- function(pch_vec) {
  u <- unique(pch_vec)
  if (length(u) == 1) u else 1
}

.draw_shape_legend <- function(d, pin, in_to_data_x, in_to_data_y,
                                shapes = field_shape_by_type(), labels = field_label_by_type()) {
  if (is.null(d$field) || !is.null(d$w1)) return(invisible(NULL))
  fields <- sort(unique(d$field))
  if (length(fields) <= 1) return(invisible(NULL))
  legend_cex <- 0.5
  x_in <- pin[1] - 0.15
  gap_in <- 0.14
  y_top_in <- pin[2] - 0.68
  for (i in seq_along(fields)) {
    y_in <- y_top_in - (i - 1) * gap_in
    pch_i <- unname(shapes[fields[i]])
    if (is.na(pch_i)) pch_i <- 16
    points(in_to_data_x(x_in), in_to_data_y(y_in), pch = pch_i, cex = 1, col = "#666666", xpd = NA)
    label_i <- unname(labels[fields[i]])
    if (is.na(label_i)) label_i <- fields[i]
    text(in_to_data_x(x_in - 0.12), in_to_data_y(y_in), label_i,
         cex = legend_cex, adj = c(1, 0.5), xpd = NA)
  }
}

# Replaces an earlier .draw_end_labels (removed 2026-10-16) that drew
# value_label/not_value_label as large text INSIDE the data panel's
# own corners, which collided with points/badges. This instead draws
# one real double-headed arrow (arrows(code=3), not a plotmath %<->%
# character -- user explicitly wants a drawn line, not a glyph, 2026-
# 10-16) per axis, spanning that axis's own FULL length (usr[1] to
# usr[2] for X, usr[3] to usr[4] for Y), placed
# in the margin between the tick labels and that axis's own title.
# Margin "line" units are converted to data coordinates via par("csi")
# (inches per line, independent of cex/canvas_scale) and the panel's
# own inches-per-data-unit ratio (same ratio panel_coord_converters()
# computes, inlined here since this runs in the margin, outside usr,
# where panel_coord_converters's own in_to_data_*(positive inches)
# convention isn't what's needed -- this wants negative/leftward
# offsets instead).
# X axis keeps value_label/not_value_label (point->target signed KLD, a
# genuinely bipolar "F" vs "not F" quantity). Y axis is fixed
# "confident"/"inconfident" wording (2026-11-05, user-specified) --
# Y is adjusted total frequency (the original node_size column, see
# draw_point_vs_target()'s own swap-block comment and the file header
# comment): a low value means the point/pair itself is rare (sparse
# evidence, so whatever its X-axis KLD reading says is less trustworthy),
# a high value means it is well-attested (frequent, so its KLD reading
# rests on more evidence) -- "inconfident" goes at the bottom (low
# adjusted frequency, near the axis origin), "confident" at the top
# (high adjusted frequency), same fixed-wording pattern dispersion
# used to have here ("uniform"/"concentrated") before dispersion moved
# to node size.
.draw_axis_direction_arrows <- function(value_label, not_value_label) {
  usr <- par("usr"); pin <- par("pin")
  x_per_in <- diff(usr[1:2]) / pin[1]
  y_per_in <- diff(usr[3:4]) / pin[2]
  line_in <- par("csi")
  mgp <- par("mgp")
  # Three separate rows now, each its own margin line, outward from the
  # plot: tick labels (occupy roughly [mgp[2], mgp[2]+1]) -> arrow (at
  # arrow_line) -> end-text labels (at label_line, a further 0.5 lines
  # out so the text doesn't sit on top of the arrowhead) -> axis title
  # (furthest out, at .draw_axes's own title_line = mgp[1] + 1.4) --
  # see .open_canvas's own comment for the matching mar increase that
  # makes room for all three new rows.
  arrow_line <- mgp[2] + 1.4
  label_line <- arrow_line + 0.5 * .text_scale()
  label_cex <- 0.8 * .text_scale()

  hide_x <- isTRUE(getOption("kld.plot.hide_x_title", FALSE))
  hide_y <- isTRUE(getOption("kld.plot.hide_y_title", FALSE))
  if (!hide_x) {
  x_arrow_y <- usr[3] - arrow_line * line_in * y_per_in
  x_label_y <- usr[3] - label_line * line_in * y_per_in
  arrows(usr[1], x_arrow_y, usr[2], x_arrow_y, code = 3, length = 0.05,
         angle = 25, col = "#888888", lwd = 0.8, xpd = NA)
  text(usr[1], x_label_y, not_value_label, adj = c(0, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
  text(usr[2], x_label_y, value_label, adj = c(1, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
  }
  if (hide_y) return(invisible())

  y_arrow_x <- usr[1] - arrow_line * line_in * x_per_in
  y_label_x <- usr[1] - label_line * line_in * x_per_in
  arrows(y_arrow_x, usr[3], y_arrow_x, usr[4], code = 3, length = 0.05,
         angle = 25, col = "#888888", lwd = 0.8, xpd = NA)
  text(y_label_x, usr[3], "inconfident", srt = 90, adj = c(0, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
  text(y_label_x, usr[4], "confident", srt = 90, adj = c(1, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
}


resolve_annotation_inputs <- function(annotate_words, annotate_words_file,
                                       highlight, highlight_file,
                                       highlight_pairs, highlight_pairs_file,
                                       annotate_pairs = NULL, annotate_pairs_file = NULL) {
  if (!is.null(annotate_words_file)) {
    raw <- readLines(annotate_words_file, encoding = "UTF-8")
    annotate_words <- unique(raw[nchar(raw) > 0])
  }
  if (!is.null(highlight_file)) {
    highlight <- read.table(highlight_file, sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  if (!is.null(highlight_pairs_file)) {
    highlight_pairs <- read.table(highlight_pairs_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  if (!is.null(annotate_pairs_file)) {
    annotate_pairs <- read.table(annotate_pairs_file, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  list(annotate_words = annotate_words, highlight = highlight, highlight_pairs = highlight_pairs,
       annotate_pairs = annotate_pairs)
}

# Named by every column either the word- or pair-level visualization
# tsv schema can carry (see load_and_merge_input's own fread() call
# below) -- fread only applies colClasses entries whose name matches
# an actual column in THIS file; names from the OTHER schema just
# aren't present in any one file, which fread warns about (benign,
# suppressWarnings()'d below) rather than erroring on.
.POINT_TABLE_COL_CLASSES <- c(
  point = "character", frequency_type = "character", kld_variant = "character",
  keyness = "numeric", frequency = "numeric", total_freq = "numeric",
  n_target = "numeric", n_not_target = "numeric", node_size = "numeric",
  dispersion = "numeric", dispersion_difference = "numeric", field = "character",
  w1 = "character", w2 = "character", field1 = "character", field2 = "character",
  freq1 = "numeric", freq2 = "numeric", dispersion1 = "numeric", dispersion2 = "numeric",
  assoc_w1_to_w2 = "numeric", assoc_w2_to_w1 = "numeric",
  keyness1 = "numeric", keyness2 = "numeric"
)

load_and_merge_input <- function(input) {
  # inputはカンマ区切りで複数パス指定可——pair-entry-entry/pair-entry-concept/
  # pair-entry-topicのような複数のenriched tsvを1枚の図にまとめて描く時に使う。
  # 列が完全に揃っていない場合(field1/field2がまだ無い古いキャッシュtsvなど)
  # でも落ちないよう.rbind_fill()で欠けている列はNA埋めする——単一ファイル
  # 指定という既存の使い方は今までどおり(splitしても要素数1、rbind_fillも
  # 1個だけならそのまま返す)。
  #
  # data.table::fread(), not read.table() -- 2026-10-17, confirmed
  # against the real pair-entry-entry-vs-anthology-shinkokin.tsv
  # (1.8M rows): read.table (even tuned with encoding=/colClasses) was
  # still ~5s on that file, too slow for Shiny's interactive adjust-
  # and-preview loop; fread is built for exactly this. data.table =
  # FALSE -- every downstream filter_by_*/subset in this file does
  # plain d[keep, ] data.frame-style row filtering, which a
  # data.table's own `[` would interpret with different (NSE)
  # semantics; converting back to a plain data.frame right at read
  # time keeps every existing call site correct unchanged.
  input_paths <- strsplit(input, ",", fixed = TRUE)[[1]]
  tables <- lapply(input_paths, function(p) {
    suppressWarnings(data.table::fread(
      p, sep = "\t", header = TRUE, quote = "", encoding = "UTF-8",
      na.strings = c("NA", ""), colClasses = .POINT_TABLE_COL_CLASSES,
      data.table = FALSE
    ))
  })
  if (length(tables) == 1) tables[[1]] else .rbind_fill(tables)
}

compute_point_shapes <- function(d, point_field) {
  if (is.null(d$w1) && !is.null(d$field)) {
    pch_vec <- unname(field_shape_by_type()[d$field])
    pch_vec[is.na(pch_vec)] <- 1
    pch_vec
  } else {
    rep(if (identical(point_field, "entry")) 1 else 0, nrow(d))
  }
}

select_badges <- function(d, annotate_words, label_top_n, neutral, annotate_nontarget_only,
                           annotate_point_exclude, badge_selection, frontier_n,
                           aggregate_locators, annotate_pairs = NULL) {
  rows <- .select_badge_rows(d, annotate_words, label_top_n, neutral,
                              annotate_nontarget_only, annotate_point_exclude,
                              badge_selection = badge_selection, frontier_n = frontier_n,
                              annotate_pairs = annotate_pairs)
  if (aggregate_locators) rows <- .aggregate_locator_groups(rows)
  rows
}

panel_coord_converters <- function() {
  usr <- par("usr")
  pin <- par("pin")
  x_per_in <- (usr[2] - usr[1]) / pin[1]
  y_per_in <- (usr[4] - usr[3]) / pin[2]
  list(
    usr = usr, pin = pin, x_per_in = x_per_in, y_per_in = y_per_in,
    data_to_in_x = function(x) (x - usr[1]) / x_per_in,
    data_to_in_y = function(y) (y - usr[3]) / y_per_in,
    in_to_data_x = function(x) usr[1] + x * x_per_in,
    in_to_data_y = function(y) usr[3] + y * y_per_in
  )
}



should_hide_unannotated_in_zoom <- function(zoom_annotated_only) {
  isTRUE(zoom_annotated_only)
}

should_hide_unannotated_for_repel <- function(repel_annotated_only, badge_layout) {
  isTRUE(repel_annotated_only) && identical(badge_layout, "repel")
}

draw_point_vs_target <- function(input, output, point_label = "point",
                                  axis_label = "target", value_label = "value",
                                  not_value_label = paste0("¬", value_label),
                                  label_top_n = 10, nofreq = FALSE,
                                  annotate_words = NULL,
                                  annotate_words_file = NULL,
                                  highlight = NULL,
                                  highlight_file = NULL,
                                  highlight_pairs = NULL,
                                  highlight_pairs_file = NULL,
                                  annotate_nontarget_only = FALSE,
                                  annotate_point_exclude = NULL,
                                  highlight_frontier = FALSE,
                                  default_label = point_label,
                                  badge_layout = "empty_region",
                                  badge_selection = "words",
                                  frontier_n = 15,
                                  canvas_size = 450,
                                  inset = FALSE,
                                  inset_position = "top",
                                  inset_configs = NULL,
                                  inset_border_col = "#FF073A",
                                  inset_pad_frac = 0.15,
                                  zoom_annotated_only = FALSE,
                                  repel_annotated_only = FALSE,
                                  edge_threshold = NULL,
                                  edge_threshold_op = "gte",
                                  dispersion_threshold = NULL,
                                  dispersion_threshold_op = "gte",
                                  node_dispersion_threshold = NULL,
                                  node_dispersion_threshold_op = "gte",
                                  frequency_threshold = NULL,
                                  frequency_threshold_op = "gte",
                                  node_frequency_threshold = NULL,
                                  node_frequency_threshold_op = "gte",
                                  node_keyness_threshold = NULL,
                                  node_keyness_threshold_op = "gte",
                                  node_signed_keyness_threshold = NULL,
                                  node_signed_keyness_threshold_op = "gte",
                                  aggregate_locators = TRUE,
                                  max_group_members = 50,
                                  network_scale = 0.01,
                                  network_node_gap = 1,
                                  point_field = "entry",
                                  annotate_pairs = NULL,
                                  annotate_pairs_file = NULL,
                                  kld_variant = NULL,
                                  frequency_type = NULL,
                                  label_lang = "ja",
                                  unselected_alpha = 0.05,
                                  main_position = NULL,
                                  inset_header_cex = 0.75,
                                  inset_label_cex = NULL,
                                  box_label_cex = 0.8,
                                  text_scale = 1,
                                  hide_x_title = FALSE,
                                  hide_y_title = FALSE,
                                  hide_x_ticks = FALSE,
                                  hide_y_ticks = FALSE) {
  # hide_x_title/hide_y_title(2026-11-14)：外部で複数枚を格子状に並べる
  # 場合(図ex-7のdraft-short.org版)、最下段以外・最左列以外の格の軸の題と
  # 方向矢印を省く。inset=TRUEの明示グリッドは格ごとに自前で上書きする。
  old_text_opt <- options(kld.plot.text_scale = text_scale,
                          kld.plot.hide_x_title = isTRUE(hide_x_title),
                          kld.plot.hide_y_title = isTRUE(hide_y_title),
                          kld.plot.hide_x_ticks = isTRUE(hide_x_ticks),
                          kld.plot.hide_y_ticks = isTRUE(hide_y_ticks))
  on.exit(options(old_text_opt), add = TRUE)
  # badge_layout="network"splices graphviz output in as svg XML text
  # (.splice_network_glyphs below) -- that only makes sense for a real
  # vector output, never png. Rather than make the caller get this
  # right and stop() otherwise, just force the extension here, before
  # OUTPUT is used anywhere below (both the inset and non-inset paths
  # read this same local variable).
  if (identical(badge_layout, "network") && !grepl("\\.svg$", output, ignore.case = TRUE)) {
    output <- sub("\\.[^.]*$", ".svg", output)
  }
  inputs <- resolve_annotation_inputs(annotate_words, annotate_words_file,
                                       highlight, highlight_file,
                                       highlight_pairs, highlight_pairs_file,
                                       annotate_pairs, annotate_pairs_file)
  annotate_words <- inputs$annotate_words
  highlight <- inputs$highlight
  highlight_pairs <- inputs$highlight_pairs
  annotate_pairs <- inputs$annotate_pairs
  d <- load_and_merge_input(input)
  # scripts/patches/visprepare.py now always merges all of its own
  # POINT_VARIANTS ("obs"/"nofreq"/"norm"/"keyassoc") into one tsv with a
  # "kld_variant" column rather than writing one tsv per variant --
  # KLD_VARIANT, if given, is this function's own single place to pick
  # one variant's rows back out (trust upstream shape: the column is
  # assumed present whenever this is non-NULL; NULL is a complete
  # no-op, so any older tsv without the column -- or a non-point
  # reference tsv -- is unaffected).
  if (!is.null(kld_variant)) d <- d[d$kld_variant == kld_variant, ]
  # frequency_type works the same way, one axis over -- visprepare.py's
  # own FREQUENCY_TYPES ("term"/"doc") are merged into the SAME tsv
  # too, unlike kld_variant this axis also changes n_target/n_not_
  # target/total_freq (and, for pairs, assoc_w1_to_w2/assoc_w2_to_w1),
  # not just keyness/frequency -- see visprepare.py's own module
  # docstring. Still just one more row filter here either way.
  if (!is.null(frequency_type)) d <- d[d$frequency_type == frequency_type, ]
  edge_threshold <- if (!is.null(edge_threshold)) as.numeric(edge_threshold) else NULL
  edge_threshold_op <- .resolve_threshold_op(edge_threshold_op)
  dispersion_threshold <- if (!is.null(dispersion_threshold)) as.numeric(dispersion_threshold) else NULL
  dispersion_threshold_op <- .resolve_threshold_op(dispersion_threshold_op)
  node_dispersion_threshold <- if (!is.null(node_dispersion_threshold)) as.numeric(node_dispersion_threshold) else NULL
  node_dispersion_threshold_op <- .resolve_threshold_op(node_dispersion_threshold_op)
  frequency_threshold <- if (!is.null(frequency_threshold)) as.numeric(frequency_threshold) else NULL
  frequency_threshold_op <- .resolve_threshold_op(frequency_threshold_op)
  node_frequency_threshold <- if (!is.null(node_frequency_threshold)) as.numeric(node_frequency_threshold) else NULL
  node_frequency_threshold_op <- .resolve_threshold_op(node_frequency_threshold_op)
  node_keyness_threshold <- if (!is.null(node_keyness_threshold)) as.numeric(node_keyness_threshold) else NULL
  node_keyness_threshold_op <- .resolve_threshold_op(node_keyness_threshold_op)
  node_signed_keyness_threshold <- if (!is.null(node_signed_keyness_threshold)) as.numeric(node_signed_keyness_threshold) else NULL
  node_signed_keyness_threshold_op <- .resolve_threshold_op(node_signed_keyness_threshold_op)
  d <- filter_by_edge_strength(d, edge_threshold, edge_threshold_op)
  if (nrow(d) == 0) {
    stop("No points left after filtering -- edge_threshold is too strict for this dataset")
  }
  # Snapshot taken before the frequency/keyness overwrite below -- the
  # annotated-points table (app.R's table tab) wants the real tsv-shaped
  # columns, not the plotting-convenience renames. .row_id survives
  # rbind()/subsetting downstream even when row names don't (rbind()
  # silently renames colliding ones via make.unique), so it is the robust
  # join key back from tie_box_rows to this snapshot, done once at the
  # very end of this function.
  d$.row_id <- seq_len(nrow(d))
  d_raw_cols <- d
  # 横軸→point2target(旧keyness、2026-10-15から変更なし)、縦軸→総頻度(調整後)
  # に恒久的に差し替え。この「調整後の総頻度」はtsvの生のfrequency列(KLD)
  # ではなくnode_size列そのもの——log2(total_freq+1)を表生成時に全表最大値で正規化
  # 済みの値で、以前はこれで点の大きさを決めていた(.size_scale自身の旧
  # コメント参照)。点の大きさは逆にdispersion(そのpoint/pair自身の分布の
  # 均等さ)に差し替える——結果として総頻度(調整後)とdispersionが担っていた
  # 表示(軸位置と点の大きさ)を入れ替えた形になる(2026-11-05、ユーザー
  # 確認済み：改めるのは縦軸だけで横軸のkeynessはそのまま、word/pair両方
  # の表に適用。最初frequency列=KLDだと誤解して軸に出したが、ユーザーが
  # 「粗頻度の調整頻度」と明確に訂正——2026-10-17のnode_size計算の議論
  # そのものを指す。ファイル先頭のコメントも参照)。全てのフィルタ(上の
  # filter_by_*)はまだ元のdispersion/total_freq/keyness/node_size列を
  # 読む必要があるので、ここ(フィルタの後、zoom_box/ranges計算の前)が
  # 差し替えの場所——下流のzoom_box解決/badge選定/ネットワークグリフ配置は
  # 全部d$frequency(x)/d$keyness(y)/d$node_size(大きさ)を素直に読むだけ
  # なので、ここで値を入れ替えておけば無改修で動く。neutral参考線/hi-lo
  # 色分け(どちらもd$keyness>neutralという符号付きKLD前提のロジック)は
  # d$keynessが今は総頻度(調整後、非負、符号なし)を保持しているので、本来の
  # 前提とは噛み合わなくなるが、dispersionだった2026-10-15時点と同じく
  # 「既知の見た目のクセとして残す、今は直さない」というユーザー確認済み
  # の方針を踏襲する。
  adjusted_frequency <- d$node_size
  d$frequency <- d$keyness
  d$keyness <- adjusted_frequency
  d$node_size <- d$dispersion

  ranges <- .compute_axis_ranges(d, nofreq)
  breaks <- .apply_nofreq_breaks(d, nofreq)
  d <- breaks$d
  neutral <- if (!is.null(breaks$y_break)) breaks$y_break$fwd(ranges$neutral) else ranges$neutral
  selection_spec <- .selection_filter_spec(
    dispersion = dispersion_threshold, dispersion_op = dispersion_threshold_op,
    node_dispersion = node_dispersion_threshold, node_dispersion_op = node_dispersion_threshold_op,
    frequency = frequency_threshold, frequency_op = frequency_threshold_op,
    node_frequency = node_frequency_threshold, node_frequency_op = node_frequency_threshold_op,
    node_keyness = node_keyness_threshold, node_keyness_op = node_keyness_threshold_op,
    node_signed_keyness = node_signed_keyness_threshold,
    node_signed_keyness_op = node_signed_keyness_threshold_op
  )
  cand_idx <- .selection_filter_idx(d, selection_spec)
  colors <- .compute_point_colors(d, highlight, highlight_pairs, annotate_point_exclude,
                                   highlight_frontier, frontier_n, frontier_candidate_idx = cand_idx,
                                   unselected_alpha = unselected_alpha)
  d$point_col <- colors$point_col
  d$point_pch <- compute_point_shapes(d, point_field)

  if (inset) {
    configs <- .normalize_inset_configs(inset_configs, inset_position, frontier_n, selection_spec)
    bundles <- .resolve_inset_bundles(d, configs, ranges, inset_pad_frac,
                                       annotate_words, label_top_n, neutral,
                                       annotate_nontarget_only, annotate_point_exclude,
                                       badge_selection, aggregate_locators, annotate_pairs,
                                       highlight_frontier = highlight_frontier, colors = colors)
    # The table-tab's "currently annotated points" output -- union of every
    # bundle's tie_box_rows (multi-config: a point annotated by ANY inset
    # counts), joined back to the pre-overwrite snapshot (see .row_id's own
    # comment above) so callers see the real tsv columns.
    annotated_ids <- unique(unlist(lapply(bundles, function(b) b$tie_box_rows$.row_id)))
    annotated_table <- d_raw_cols[d_raw_cols$.row_id %in% annotated_ids, , drop = FALSE]
    annotated_table$.row_id <- NULL
    annotated_table <- .format_display_columns(annotated_table)
    if (isTRUE(highlight_frontier) && !is.null(colors$frontier_match_idx)) {
      union_idx <- unique(unlist(lapply(bundles, function(b) b$frontier_d_idx)))
      global_idx <- which(colors$frontier_match_idx)
      stale_idx <- global_idx[!global_idx %in% union_idx]
      if (length(stale_idx) > 0) {
        gp <- gradient_color_palette()
        d$point_col[stale_idx] <- adjustcolor(
          .linear_ramp(colors$dispersion[stale_idx], colors$signed_range, gp["lo"], gp["hi"]),
          alpha.f = colors$unselected_alpha)
      }
      if (length(union_idx) > 0) {
        d$point_col[union_idx] <- .abs_ramp(colors$dispersion[union_idx], colors$abs_range,
                                             colors$frontier_pal["light"], colors$frontier_pal["dark"])
      }
    }
    .render_with_insets(d, bundles, colors, output, ranges, breaks, neutral, nofreq,
                         axis_label, point_label, value_label, not_value_label,
                         default_label, badge_selection, badge_layout,
                         canvas_size, inset_border_col,
                         max_group_members = max_group_members,
                         network_scale = network_scale,
                         network_node_gap = network_node_gap,
                         zoom_annotated_only = zoom_annotated_only,
                         repel_annotated_only = repel_annotated_only,
                         kld_variant = kld_variant, label_lang = label_lang,
                         main_position = main_position, inset_header_cex = inset_header_cex,
                         inset_label_cex = inset_label_cex, box_label_cex = box_label_cex)
    cat(output, ":", nrow(d), point_label, "points\n")
    return(invisible(annotated_table))
  }

  d_cand <- d[cand_idx, ]
  badge_rows <- select_badges(d_cand, annotate_words, label_top_n, neutral,
                               annotate_nontarget_only, annotate_point_exclude,
                               badge_selection, frontier_n, aggregate_locators,
                               annotate_pairs = annotate_pairs)
  # .draw_x_tie_boxes自身のdocstring参照：aggregate_locators=TRUEの時、
  # badge_rowsは既に(n_target, n_not_target, kind)で1行にまとめられた後
  # (.aggregate_locator_groups)なので、同じキーでもう一度グループ化しよ
  # うとしても各グループが必ず1行になり、箱が一つも描けなくなる
  # (2026-10-17にユーザー報告で発覚)。aggregate_locators=FALSEで同じ
  # select_badgesをもう一度呼び、集約前の個々の点(「annotateされた点」
  # ではあるが、まだ合算されていない実際の座標を持つ行)を別に用意する。
  tie_box_rows <- select_badges(d_cand, annotate_words, label_top_n, neutral,
                                 annotate_nontarget_only, annotate_point_exclude,
                                 badge_selection, frontier_n,
                                 aggregate_locators = FALSE,
                                 annotate_pairs = annotate_pairs)
  # The table-tab's "currently annotated points" output -- same row set as
  # tie_box_rows, joined back to the pre-overwrite snapshot (see .row_id's
  # own comment above) so callers see the real tsv columns.
  annotated_table <- d_raw_cols[d_raw_cols$.row_id %in% tie_box_rows$.row_id, , drop = FALSE]
  annotated_table$.row_id <- NULL
  annotated_table <- .format_display_columns(annotated_table)

  .open_canvas(output, canvas_size)
  on.exit(dev.off(), add = TRUE)
  .draw_axes(d, ranges, nofreq, axis_label, point_label, breaks$x_break, breaks$y_break,
             kld_variant = kld_variant, label_lang = label_lang)

  d_drawn <- if (should_hide_unannotated_for_repel(repel_annotated_only, badge_layout)) {
    filter_to_annotated_points(d, badge_rows)
  } else d
  points(d_drawn$frequency, d_drawn$keyness, pch = d_drawn$point_pch,
         cex = .size_scale(d_drawn$node_size), col = d_drawn$point_col)

  pc <- panel_coord_converters()
  badge_cex <- 0.55 * .text_scale()
  fam <- getOption("kld.plot.font_family", "Hiragino Mincho ProN")
  label_gap_in <- strheight("あ", cex = badge_cex, family = fam) * 0.5

  tie_boxes <- .compute_x_tie_boxes(
    tie_box_rows, point_radius_in = .cex_to_radius_in(.size_scale(tie_box_rows$node_size))
  )
  .draw_x_tie_boxes(tie_boxes)
  badge_colors <- .badge_color_from_point(badge_rows)
  .draw_badges(badge_rows, fill = badge_colors$fill, border = badge_colors$border_col,
               text_col = badge_colors$text_col, word_col = badge_colors$word_col,
               badge_layout = badge_layout, d_all = d_drawn,
               value_label = value_label, not_value_label = not_value_label,
               badge_cex = badge_cex, fam = fam, label_gap_in = label_gap_in,
               x_per_in = pc$x_per_in, y_per_in = pc$y_per_in,
               data_to_in_x = pc$data_to_in_x, data_to_in_y = pc$data_to_in_y,
               in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y,
               pin = pc$pin, neutral = neutral, max_group_members = max_group_members,
               tie_boxes = tie_boxes)

  .draw_axis_direction_arrows(value_label, not_value_label)
  .draw_size_legend(d$node_size, pc$pin, pc$in_to_data_x, pc$in_to_data_y,
                     pch = .representative_pch(d$point_pch), label_lang = label_lang)
  .draw_shape_legend(d, pc$pin, pc$in_to_data_x, pc$in_to_data_y)
  .draw_legend(colors$type_ramps, colors$signed_range, colors$legend_cols, default_label, pc$pin,
               pc$in_to_data_x, pc$in_to_data_y, colors$default_palette, label_lang = label_lang)

  cat(output, ":", nrow(d), point_label, "points\n")
  invisible(annotated_table)
}


.compute_zoom_box <- function(badge_rows, main_ranges, pad_frac, min_span_frac = 0.05) {
  x_rng <- range(badge_rows$frequency); y_rng <- range(badge_rows$keyness)
  xpad <- max(diff(x_rng) * pad_frac, diff(main_ranges$xr) * min_span_frac)
  ypad <- max(diff(y_rng) * pad_frac, diff(main_ranges$yr) * min_span_frac)
  list(x1 = x_rng[1] - xpad, x2 = x_rng[2] + xpad,
       y1 = y_rng[1] - ypad, y2 = y_rng[2] + ypad)
}

.INSET_PALETTE <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7")

.draw_source_box <- function(box, border_col, label = NULL, lty = "1343", label_cex = 0.8,
                             label_x = box$x1) {
  rect(box$x1, box$y1, box$x2, box$y2, border = border_col, lty = lty,
       lwd = 1.2, xpd = NA)
  if (!is.null(label) && nzchar(label)) {
    text(label_x, box$y2, labels = label, adj = c(1.2, -0.3), xpd = NA,
         col = border_col, font = 2, cex = label_cex)
  }
}

# 2026-10-17、ユーザー確認済み：X軸(point→target KLD)の値が完全一致する
# tieでグループ化していた最初の版は、低頻度語の大半を取りこぼしていた
# ことが判明した -- pair_cross_table自身のdocstring通り、population
# (分母側)はMEASURED自身が実際に現れる行ごとに計算される(scope全体で
# 共有される定数ではない)ため、n_target/n_not_targetが同じ組み合わせ
# でも実際に現れる行が違えばKLDはわずかに(だが非ゼロに)ズレる--図の
# 上では同じ場所に密集して見えても、浮動小数点としては別の値になる。
# なのでグループ化の基準をX(keyness)の完全一致から、ROWSが実際に持つ
# n_target/n_not_target(target/neg_target側のmatch数、こちらは本当に
# 同じ値になる)の一致に変更。X・Y両方ともその組の値の範囲(min~max)を
# 取った外接矩形を描く(前版は常にX固定の細い縦長の箱だったが、今は同じ
# グループ内でもXが微妙にばらつくため、Xにも幅を持たせる必要がある)。
# MIN_GROUP_SIZE=2(2026-10-17、ユーザー確認済み：一度3に引き上げたが、
# ちょうど2点共有している組も箱で強調したいとのことで2に戻した -- 2点
# 以上共有している組は全て箱で囲む)。INSET_PAD_FRACで外接矩形の外側に
# 少し余白を持たせる(点のマーカーぎりぎりに箱の辺が重ならないように)。
# badge_selection問わず(前版はbadge_selection="zoom_box"限定だったが、
# この分析自体はbadge選択方式と無関係なので2026-10-17に限定を外した)
# 常に実行する。
#
# ROWSには実際にannotateされた点のみを渡す(d_drawn/dd--描画されている
# 点全部--を試したところ、annotateされていない点まで箱に入るのはノイズ
# が多すぎるとユーザーからフィードバックがあり、2026-10-17に絞った)が、
# aggregate_locators=TRUEの時のbadge_rowsをそのまま渡してはいけない --
# .aggregate_locator_groups()が既に同じ(n_target, n_not_target, kind)
# キーで1行にまとめた後なので、ここでもう一度同じキーでグループ化しても
# 各グループが必ず1行になり、箱が一つも描けない(2026-10-17ユーザー報告
# で発覚："框完全消失了")。呼び出し元(draw_point_vs_target)は
# aggregate_locators=FALSEでselect_badgesをもう一度呼んだ
# tie_box_rows(「annotateされた点」ではあるが集約されていない、各点
# 自身の実座標を保ったままの行)を渡すこと -- badge_rows自体(集約後、
# badge表示用)は変更しない。
#
# .compute_x_tie_boxes()はrect()を直接描かず、{key: {x1,y1,x2,y2}}
# (dataの座標系のまま)を返すだけ -- .draw_badges()自身の連結線も、この
# 点のbadgeがどこかの箱に属している時は箱の辺から出発させたい
# (2026-10-17、ユーザー確認済み：「線がもし框なら、框のedgeから出発」)
# ので、計算結果を.draw_x_tie_boxes()(実際に描く)と.draw_badges()
# (連結線の出発点を決める)の両方で共有できるよう、計算と描画を分離した。
# POINT_RADIUS_IN(ROWSと同じ長さ、各点を実際に描いた時のpch記号の半径を
# インチで)はオプション -- 渡されなければ0扱い(旧来通り、usrのごく小さい
# 一定割合だけの余白)。渡す理由(2026-10-17ユーザー報告：拡大図で枠が点の
# 拡大に追従せず、点が枠からはみ出す)：INSET_PAD_FRACはusrの範囲に対する
# 固定割合なので、インチ換算すると「パネル幅の一定%」という、ズーム
# レベルに関係なく常に同じ絶対サイズになる(usrが狭くなるほど1データ単位
# あたりのインチ数が増えるので、data単位の余白をインチに戻すと両者が
# 打ち消し合う)。一方render_zoomed_scatter_panel側の点のcexはzoom_factor
# (最大4倍)で別途拡大されるため、ズームパネルでは記号の実際の半径が
# パネル幅に対して不釣り合いに大きくなり、この固定インチ幅の余白では
# 包みきれなくなる。各グループの最大記号半径(呼び出し元がそのパネルで
# 実際に使ったcex式と同じ式で算出し、インチのまま渡す)をデータ単位に
# 換算して、INSET_PAD_FRAC由来の余白とのmaxを取ることで、記号がどれだけ
# 拡大されていても箱が必ず記号を覆うようにする。
.compute_x_tie_boxes <- function(rows, point_radius_in = NULL, inset_pad_frac = 0.01,
                                  min_group_size = 2) {
  if (nrow(rows) < min_group_size) return(list())
  usr <- par("usr"); pin <- par("pin")
  xpad <- diff(usr[1:2]) * inset_pad_frac
  ypad <- diff(usr[3:4]) * inset_pad_frac
  x_per_in <- diff(usr[1:2]) / pin[1]
  y_per_in <- diff(usr[3:4]) / pin[2]
  if (is.null(point_radius_in)) point_radius_in <- rep(0, nrow(rows))
  key <- paste(rows$n_target, rows$n_not_target, sep = "/")
  groups <- split(seq_len(nrow(rows)), key)
  boxes <- list()
  for (k in names(groups)) {
    idx <- groups[[k]]
    if (length(idx) < min_group_size) next
    x_rng <- range(rows$frequency[idx])
    y_rng <- range(rows$keyness[idx])
    if (diff(x_rng) <= 0 && diff(y_rng) <= 0) next
    max_radius_in <- max(point_radius_in[idx], 0)
    pad_x <- max(xpad, max_radius_in * x_per_in)
    pad_y <- max(ypad, max_radius_in * y_per_in)
    boxes[[k]] <- list(x1 = x_rng[1] - pad_x, y1 = y_rng[1] - pad_y,
                        x2 = x_rng[2] + pad_x, y2 = y_rng[2] + pad_y)
  }
  boxes
}

# cex(文字/記号サイズの倍率)から、その記号(pch 0:25系)がおおよそ占める
# 半径をインチで概算する -- Rのpch記号は1文字分の高さ(par("cin")[2])を
# 基準にcex倍されるので、直径をcex*cin[2]、半径をその半分とする近似。
# ピクセル単位で厳密に一致させる必要はない(.compute_x_tie_boxes側は
# これにinset_pad_fracの余白もさらに足すので、多少大きめに見積もる分には
# 問題なく、小さすぎる方だけが実害になる)。
.cex_to_radius_in <- function(cex) cex * par("cin")[2] * 0.5

# 単純な箱のみ描く(badge=文字ラベル付きではない)、グループを選別する
# 閾値/上限は.compute_x_tie_boxes側のMIN_GROUP_SIZE以外に設けない
# (ユーザー確認済み)。灰色の破線(2026-10-17、赤の実線から変更、
# ユーザー確認済み)。
.draw_x_tie_boxes <- function(boxes, border_col = "#999999") {
  for (box in boxes) {
    rect(box$x1, box$y1, box$x2, box$y2, border = border_col, lty = 2, lwd = 1, xpd = NA)
  }
  invisible(NULL)
}

.draw_projection_lines <- function(src_ndc_x, src_ndc_y, dest_ndc_x, dest_ndc_y, col) {
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  corners <- expand.grid(xi = 1:2, yi = 1:2)
  for (k in seq_len(nrow(corners))) {
    xi <- corners$xi[k]; yi <- corners$yi[k]
    segments(src_ndc_x[xi], src_ndc_y[yi], dest_ndc_x[xi], dest_ndc_y[yi],
             col = col, lty = "1343", lwd = 0.8, xpd = NA)
  }
}

render_inset_source_panel <- function(d, ranges, nofreq, axis_label, point_label, breaks,
                                       zboxes, box_cols, value_label, not_value_label,
                                       colors, default_label, kld_variant = NULL,
                                       box_labels, box_ltys, label_lang = "ja",
                                       box_label_cex = 0.8) {
  .draw_axes(d, ranges, nofreq, axis_label, point_label, breaks$x_break, breaks$y_break,
             kld_variant = kld_variant, label_lang = label_lang)
  points(d$frequency, d$keyness, pch = d$point_pch,
         cex = .size_scale(d$node_size), col = d$point_col)
  src_ndc_x <- vector("list", length(zboxes))
  src_ndc_y <- vector("list", length(zboxes))
  # 字母の重なり回避(2026-11-14)：ほぼ同じ框角に複数の字母が来る場合
  # (例：閾値の近いD/Eの框がほぼ一致)、先に置いた字母と重なる分だけ
  # 左へずらす。重ならない場合はlabel_x=x1のままで従来と同じ位置。
  placed <- list()
  for (i in seq_along(zboxes)) {
    zb <- zboxes[[i]]
    lab <- box_labels[i]
    label_x <- zb$x1
    if (!is.null(lab) && nzchar(lab)) {
      lw <- strwidth(lab, cex = box_label_cex, font = 2) * 1.2
      lh <- strheight(lab, cex = box_label_cex, font = 2) * 1.2
      repeat {
        hit <- vapply(placed, function(p) abs(p[1] - label_x) < lw && abs(p[2] - zb$y2) < lh,
                      logical(1))
        if (!any(hit)) break
        label_x <- label_x - lw
      }
      placed[[length(placed) + 1]] <- c(label_x, zb$y2)
    }
    .draw_source_box(zb, box_cols[i], label = lab, lty = box_ltys[i],
                     label_cex = box_label_cex, label_x = label_x)
    src_ndc_x[[i]] <- grconvertX(c(zboxes[[i]]$x1, zboxes[[i]]$x2), from = "user", to = "ndc")
    src_ndc_y[[i]] <- grconvertY(c(zboxes[[i]]$y1, zboxes[[i]]$y2), from = "user", to = "ndc")
  }
  pc <- panel_coord_converters()
  .draw_axis_direction_arrows(value_label, not_value_label)
  .draw_size_legend(d$node_size, pc$pin, pc$in_to_data_x, pc$in_to_data_y,
                     pch = .representative_pch(d$point_pch), label_lang = label_lang)
  .draw_shape_legend(d, pc$pin, pc$in_to_data_x, pc$in_to_data_y)
  .draw_legend(colors$type_ramps, colors$signed_range, colors$legend_cols, default_label, pc$pin,
               pc$in_to_data_x, pc$in_to_data_y, colors$default_palette, label_lang = label_lang)
  list(src_ndc_x = src_ndc_x, src_ndc_y = src_ndc_y)
}

render_network_glyph_panel <- function(badge_rows, inset_border_col, value_label, not_value_label,
                                        network_scale, network_node_gap) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  box(col = inset_border_col, lwd = 1.5)
  pc <- panel_coord_converters()
  dest_ndc_x <- grconvertX(pc$usr[1:2], from = "user", to = "ndc")
  dest_ndc_y <- grconvertY(pc$usr[3:4], from = "user", to = "ndc")

  badge_cex_net <- 0.55 * .text_scale()
  fam_net <- getOption("kld.plot.font_family", "Hiragino Mincho ProN")
  count_label_net <- sprintf("%s:%d / %s:%d", value_label, badge_rows$n_target[1],
                              not_value_label, badge_rows$n_not_target[1])
  count_w_in <- strwidth(count_label_net, cex = badge_cex_net, family = fam_net, units = "inches")
  count_h_in <- strheight(count_label_net, cex = badge_cex_net, family = fam_net, units = "inches")
  net_box_w_in <- count_w_in + count_w_in * 0.35
  net_box_h_in <- count_h_in + count_h_in * 1.2
  net_pad_in <- 0.08
  net_colors <- .badge_color_from_point(badge_rows[1, , drop = FALSE])
  rect(pc$in_to_data_x(net_pad_in), pc$in_to_data_y(pc$pin[2] - net_pad_in - net_box_h_in),
       pc$in_to_data_x(net_pad_in + net_box_w_in), pc$in_to_data_y(pc$pin[2] - net_pad_in),
       col = net_colors$fill, border = net_colors$border_col, xpd = NA)
  text(pc$in_to_data_x(net_pad_in + net_box_w_in / 2),
       pc$in_to_data_y(pc$pin[2] - net_pad_in - net_box_h_in / 2),
       labels = count_label_net, cex = badge_cex_net, family = fam_net,
       col = net_colors$text_col, xpd = NA)

  members <- if (!is.null(badge_rows$group_size)) {
    badge_rows$members[[1]]
  } else {
    net_member_cols <- c("w1", "w2", "freq1", "freq2", "assoc_w1_to_w2", "assoc_w2_to_w1")
    if (!is.null(badge_rows$field1)) net_member_cols <- c(net_member_cols, "field1", "field2")
    badge_rows[, net_member_cols]
  }
  svg_txt <- .render_pair_network_svg(members, network_scale = network_scale,
                                        network_node_gap = network_node_gap)
  dims <- .parse_svg_dims_pt(svg_txt)
  pad_in <- 0.15
  avail_w_in <- max(pc$pin[1] - 2 * pad_in, 0.1)
  avail_h_in <- max(pc$pin[2] - 2 * pad_in, 0.1)
  fit_scale <- min(avail_w_in * 72 / dims$w, avail_h_in * 72 / dims$h)
  w_in <- dims$w * fit_scale / 72
  h_in <- dims$h * fit_scale / 72
  pending <- list(list(svg_txt = svg_txt, gv_w_pt = dims$w, gv_h_pt = dims$h,
                        x1_in = (pc$pin[1] - w_in) / 2, y1_in = (pc$pin[2] - h_in) / 2,
                        w_in = w_in, h_in = h_in))
  list(pending = pending, dest_ndc_x = dest_ndc_x, dest_ndc_y = dest_ndc_y,
       in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y, pin = pc$pin)
}

render_zoomed_scatter_panel <- function(d, badge_rows, tie_box_rows, zbox, ranges, nofreq, axis_label, point_label,
                                         inset_border_col, badge_selection, badge_layout,
                                         value_label, not_value_label, neutral, max_group_members,
                                         zoom_annotated_only = FALSE, repel_annotated_only = FALSE,
                                         kld_variant = NULL, label_lang = "ja") {
  box_ranges <- list(xr = c(zbox$x1, zbox$x2), yr = c(zbox$y1, zbox$y2),
                      xpad = 0, ypad = 0, neutral = ranges$neutral)
  .draw_axes(d, box_ranges, nofreq, axis_label, point_label, NULL, NULL, axis_col = inset_border_col,
             kld_variant = kld_variant, label_lang = label_lang, guides = FALSE)
  dd <- d[d$frequency >= zbox$x1 & d$frequency <= zbox$x2 &
            d$keyness >= zbox$y1 & d$keyness <= zbox$y2, ]
  if (should_hide_unannotated_in_zoom(zoom_annotated_only) ||
        should_hide_unannotated_for_repel(repel_annotated_only, badge_layout)) {
    dd <- filter_to_annotated_points(dd, badge_rows)
  }
  zoom_x <- (diff(ranges$xr) + 2 * ranges$xpad) / diff(c(zbox$x1, zbox$x2))
  zoom_y <- (diff(ranges$yr) + 2 * ranges$ypad) / diff(c(zbox$y1, zbox$y2))
  zoom_factor <- .capped_zoom_factor(sqrt(zoom_x * zoom_y))
  points(dd$frequency, dd$keyness, pch = dd$point_pch,
         cex = zoom_factor * .size_scale(dd$node_size), col = dd$point_col)
  box(col = inset_border_col, lwd = 1.5)

  pc <- panel_coord_converters()
  badge_cex <- 0.55 * .text_scale()
  fam <- getOption("kld.plot.font_family", "Hiragino Mincho ProN")
  label_gap_in <- strheight("あ", cex = badge_cex, family = fam) * 0.5
  dest_ndc_x <- grconvertX(pc$usr[1:2], from = "user", to = "ndc")
  dest_ndc_y <- grconvertY(pc$usr[3:4], from = "user", to = "ndc")

  tie_boxes <- .compute_x_tie_boxes(
    tie_box_rows,
    point_radius_in = .cex_to_radius_in(zoom_factor * .size_scale(tie_box_rows$node_size))
  )
  .draw_x_tie_boxes(tie_boxes)
  badge_colors <- .badge_color_from_point(badge_rows)
  .draw_badges(badge_rows, fill = badge_colors$fill, border = badge_colors$border_col,
               text_col = badge_colors$text_col, word_col = badge_colors$word_col,
               badge_layout = badge_layout, d_all = dd,
               value_label = value_label, not_value_label = not_value_label,
               badge_cex = badge_cex, fam = fam, label_gap_in = label_gap_in,
               x_per_in = pc$x_per_in, y_per_in = pc$y_per_in,
               data_to_in_x = pc$data_to_in_x, data_to_in_y = pc$data_to_in_y,
               in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y,
               pin = pc$pin, neutral = neutral, max_group_members = max_group_members,
               tie_boxes = tie_boxes)
  .draw_axis_direction_arrows(value_label, not_value_label)
  list(pending = list(), dest_ndc_x = dest_ndc_x, dest_ndc_y = dest_ndc_y,
       in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y, pin = pc$pin)
}

.INSET_GRID_CELLS <- list(
  top_left = c(col = 1, row = 3), top = c(col = 2, row = 3), top_right = c(col = 3, row = 3),
  left     = c(col = 1, row = 2),                            right     = c(col = 3, row = 2),
  bottom_left = c(col = 1, row = 1), bottom = c(col = 2, row = 1), bottom_right = c(col = 3, row = 1)
)

.GRID_CELL_RE <- "^r([0-9]+)c([0-9]+)$"
.GRID_MAR_FULL <- 5.0
.GRID_MAR_TICKS <- 1.5
.GRID_MAR_RIGHT <- 1.0 # 左右の間隔はやや広め、上下は詰める(2026-11-14)
.explicit_grid_layout <- function(positions, main_position, top_lines = 2.6) {
  if (is.null(main_position)) main_position <- "r1c1"
  all_pos <- c(main_position, positions)
  bad <- all_pos[!grepl(.GRID_CELL_RE, all_pos)]
  if (length(bad) > 0) {
    stop("显式网格(main_position/rRcC)时所有inset_position都必须是r<行>c<列>写法，收到: ",
         paste(bad, collapse = ", "))
  }
  if (length(all_pos) != length(unique(all_pos))) {
    stop("主图/inset的网格位置有重复: ", paste(all_pos, collapse = ", "))
  }
  rc <- regmatches(all_pos, regexec(.GRID_CELL_RE, all_pos))
  rr <- as.integer(vapply(rc, `[`, "", 2)); cc <- as.integer(vapply(rc, `[`, "", 3))
  nr <- max(rr); nc <- max(cc)
  cell_lines <- 300 / 72 / 0.2
  square <- cell_lines - 6
  left <- ifelse(seq_len(nc) == 1, .GRID_MAR_FULL, .GRID_MAR_TICKS)
  bottom <- ifelse(seq_len(nr) == nr, .GRID_MAR_FULL, .GRID_MAR_TICKS)
  col_w <- square + left + .GRID_MAR_RIGHT
  row_h <- square + bottom + top_lines
  W <- sum(col_w); H <- sum(row_h)
  col_x0 <- c(0, cumsum(col_w))[seq_len(nc)]
  row_top <- c(0, cumsum(row_h))[seq_len(nr)]
  cell <- function(r, c) {
    # 浮動小数の誤差で0/1をわずかに外れるとpar(fig=)がエラーになるので丸める
    fig <- c(col_x0[c] / W, (col_x0[c] + col_w[c]) / W,
             (H - row_top[r] - row_h[r]) / H, (H - row_top[r]) / H)
    list(fig = pmin(pmax(round(fig, 10), 0), 1),
         mar = c(bottom[r], left[c], top_lines, .GRID_MAR_RIGHT),
         hide_x = r != nr, hide_y = c != 1)
  }
  cells <- lapply(seq_along(all_pos), function(i) cell(rr[i], cc[i]))
  list(main = cells[[1]]$fig, insets = lapply(cells[-1], `[[`, "fig"),
       main_cell = cells[[1]], inset_cells = cells[-1],
       width_multiplier = W / cell_lines, height_multiplier = H / cell_lines)
}

.inset_grid_layout <- function(positions, main_position = NULL, top_lines = 2.6) {
  if (!is.null(main_position) || any(grepl(.GRID_CELL_RE, positions))) {
    return(.explicit_grid_layout(positions, main_position, top_lines))
  }
  if (length(positions) != length(unique(positions))) {
    stop("inset_configs里有重复的inset_position: ", paste(positions, collapse = ", "))
  }
  unknown <- setdiff(positions, names(.INSET_GRID_CELLS))
  if (length(unknown) > 0) {
    stop("inset_position必须是以下之一: ", paste(names(.INSET_GRID_CELLS), collapse = "/"),
         "，收到未知值: ", paste(unknown, collapse = ", "))
  }
  cells <- .INSET_GRID_CELLS[positions]
  used_cols <- sort(unique(c(2, vapply(cells, `[[`, numeric(1), "col"))))
  used_rows <- sort(unique(c(2, vapply(cells, `[[`, numeric(1), "row"))))
  col_idx <- setNames(seq_along(used_cols), used_cols)
  row_idx <- setNames(seq_along(used_rows), used_rows)
  nc <- length(used_cols); nr <- length(used_rows)
  cell_fig <- function(col, row) {
    c1 <- col_idx[[as.character(col)]]; r1 <- row_idx[[as.character(row)]]
    c((c1 - 1) / nc, c1 / nc, (r1 - 1) / nr, r1 / nr)
  }
  main <- cell_fig(2, 2)
  insets <- lapply(cells, function(cell) cell_fig(cell["col"], cell["row"]))
  list(main = main, insets = insets, width_multiplier = nc, height_multiplier = nr)
}

.inset_fig_regions <- function(inset_position) {
  lay <- .inset_grid_layout(inset_position)
  list(main = lay$main, inset = lay$insets[[1]],
       width_multiplier = lay$width_multiplier, height_multiplier = lay$height_multiplier)
}

.normalize_inset_configs <- function(configs, inset_position, frontier_n, selection_spec) {
  if (is.null(configs) || length(configs) == 0) {
    return(list(list(label = "", inset_position = inset_position, frontier_n = frontier_n,
                      filters = selection_spec)))
  }
  n <- length(configs)
  labels <- if (n == 1) "" else letters[seq_len(n)]
  lapply(seq_len(n), function(i) {
    cfg <- configs[[i]]
    pos <- cfg$inset_position; if (is.null(pos)) pos <- inset_position
    fn <- cfg$frontier_n; if (is.null(fn)) fn <- frontier_n
    filt <- cfg$filters; if (is.null(filt)) filt <- selection_spec
    list(label = labels[i], inset_position = pos, frontier_n = fn, filters = filt)
  })
}

.resolve_inset_bundles <- function(d, configs, ranges, inset_pad_frac,
                                    annotate_words, label_top_n, neutral,
                                    annotate_nontarget_only, annotate_point_exclude,
                                    badge_selection, aggregate_locators, annotate_pairs,
                                    highlight_frontier = FALSE, colors = NULL) {
  lapply(configs, function(cfg) {
    cand_idx <- .selection_filter_idx(d, cfg$filters)
    if (length(cand_idx) == 0) {
      stop(sprintf("inset%s: 这个config的filter把候选池筛空了，没有点可以画放大图",
                    if (nzchar(cfg$label)) paste0(" ", cfg$label) else ""))
    }
    d_eligible <- d[cand_idx, ]
    frontier_d_idx <- integer(0)
    if (isTRUE(highlight_frontier) && !is.null(colors) && !is.null(colors$frontier_match_idx)) {
      was_global_frontier <- colors$frontier_match_idx[cand_idx]
      if (any(was_global_frontier)) {
        reset_idx <- cand_idx[was_global_frontier]
        gp <- gradient_color_palette()
        base_col <- adjustcolor(
          .linear_ramp(colors$dispersion[reset_idx], colors$signed_range, gp["lo"], gp["hi"]),
          alpha.f = colors$unselected_alpha)
        d_eligible$point_col[was_global_frontier] <- base_col
      }
      local_match <- .frontier_match_idx(d, cfg$frontier_n, cand_idx, annotate_point_exclude)
      local_match_in_eligible <- local_match[cand_idx]
      if (any(local_match_in_eligible)) {
        local_idx <- cand_idx[local_match_in_eligible]
        frontier_d_idx <- local_idx
        d_eligible$point_col[local_match_in_eligible] <- .abs_ramp(
          colors$dispersion[local_idx], colors$abs_range,
          colors$frontier_pal["light"], colors$frontier_pal["dark"])
      }
    }
    badge_rows <- select_badges(d_eligible, annotate_words, label_top_n, neutral,
                                 annotate_nontarget_only, annotate_point_exclude,
                                 badge_selection, cfg$frontier_n, aggregate_locators,
                                 annotate_pairs = annotate_pairs)
    tie_box_rows <- select_badges(d_eligible, annotate_words, label_top_n, neutral,
                                   annotate_nontarget_only, annotate_point_exclude,
                                   badge_selection, cfg$frontier_n, aggregate_locators = FALSE,
                                   annotate_pairs = annotate_pairs)
    zbox <- .compute_zoom_box(tie_box_rows, ranges, inset_pad_frac)
    list(label = cfg$label, inset_position = cfg$inset_position, frontier_n = cfg$frontier_n,
         desc = .inset_filter_desc(cfg$filters), d_eligible = d_eligible,
         badge_rows = badge_rows, tie_box_rows = tie_box_rows, zbox = zbox,
         frontier_d_idx = frontier_d_idx)
  })
}

.is_single_locator <- function(badge_rows) {
  is_pair <- !is.null(badge_rows$w1)
  if (!is_pair || nrow(badge_rows) == 0) return(FALSE)
  kind <- if (!is.null(badge_rows$field1)) paste(badge_rows$field1, badge_rows$field2, sep = "|") else ""
  key <- paste(badge_rows$n_target, badge_rows$n_not_target, kind, sep = "\t")
  length(unique(key)) == 1
}

.draw_inset_header <- function(label, desc, col, cex = 0.75, label_cex = NULL,
                                fam = getOption("kld.plot.font_family", "Hiragino Mincho ProN")) {
  if (is.null(label_cex) || !nzchar(label)) {
    txt <- if (nzchar(label)) sprintf("%s: %s", label, desc) else desc
    mtext(txt, side = 3, line = 0.3, adj = 0, cex = cex, col = col, family = fam)
    return(invisible())
  }
  usr <- par("usr")
  mtext(label, side = 3, line = 0.3, adj = 0, at = usr[1], cex = label_cex, font = 2,
        col = col, family = fam)
  lab_w <- strwidth(paste0(label, " "), units = "user", cex = label_cex, font = 2, family = fam)
  mtext(desc, side = 3, line = 0.3, adj = 0, at = usr[1] + lab_w, cex = cex,
        col = col, family = fam)
}

.render_with_insets <- function(d, bundles, colors, output, ranges, breaks, neutral, nofreq,
                                 axis_label, point_label, value_label, not_value_label,
                                 default_label, badge_selection, badge_layout,
                                 canvas_size, inset_border_col,
                                 max_group_members = 50, network_scale = 0.01,
                                 network_node_gap = 1,
                                 zoom_annotated_only = FALSE, repel_annotated_only = FALSE,
                                 kld_variant = NULL, label_lang = "ja",
                                 main_position = NULL, inset_header_cex = 0.75,
                                 inset_label_cex = NULL, box_label_cex = 0.8) {
  multi <- length(bundles) > 1
  single_locator_flags <- vapply(bundles, function(b) .is_single_locator(b$badge_rows), logical(1))

  positions <- vapply(bundles, `[[`, "", "inset_position")
  top_lines <- (0.2 + max(inset_header_cex, inset_label_cex)) / (canvas_size / 300)
  lay <- .inset_grid_layout(positions, main_position, top_lines)
  explicit <- !is.null(lay$main_cell)
  if (explicit) {
    old_opts <- options(kld.plot.hide_x_title = FALSE, kld.plot.hide_y_title = FALSE)
    on.exit(options(old_opts), add = TRUE)
  }
  canvas <- .open_canvas(output, canvas_size, width_multiplier = lay$width_multiplier,
                          height_multiplier = lay$height_multiplier)
  splice_ok <- FALSE
  on.exit(if (!splice_ok) dev.off(), add = TRUE)

  par(fig = lay$main)
  if (explicit) {
    par(mar = lay$main_cell$mar * canvas$canvas_scale)
    options(kld.plot.hide_x_title = lay$main_cell$hide_x,
            kld.plot.hide_y_title = lay$main_cell$hide_y)
  }
  zboxes <- lapply(bundles, `[[`, "zbox")
  box_labels <- vapply(bundles, `[[`, "", "label")
  lty_cycle <- c("1343", "2262", "4212", "73", "13")
  box_ltys <- lty_cycle[(seq_along(zboxes) - 1) %% length(lty_cycle) + 1]
  # Multi-config: one color per box, shared by the source box, its letter and
  # the matching inset's border/axes/header. Color plus the lty cycle keeps
  # boxes distinguishable for color-blind readers and in grayscale print.
  # Single config keeps inset_border_col.
  box_cols <- if (multi) {
    .INSET_PALETTE[(seq_along(zboxes) - 1) %% length(.INSET_PALETTE) + 1]
  } else {
    rep(inset_border_col, length(zboxes))
  }
  source_panel <- render_inset_source_panel(d, ranges, nofreq, axis_label, point_label, breaks,
                                             zboxes, box_cols, value_label, not_value_label,
                                             colors, default_label, kld_variant = kld_variant,
                                             box_labels = box_labels, box_ltys = box_ltys,
                                             label_lang = label_lang,
                                             box_label_cex = box_label_cex)

  all_pending <- list()
  all_marks_in <- list()
  single_dest_ndc_x <- NULL
  single_dest_ndc_y <- NULL
  for (i in seq_along(bundles)) {
    b <- bundles[[i]]
    par(fig = lay$insets[[i]], new = TRUE)
    if (multi) {
      inset_mar <- if (explicit) {
        options(kld.plot.hide_x_title = lay$inset_cells[[i]]$hide_x,
                kld.plot.hide_y_title = lay$inset_cells[[i]]$hide_y)
        lay$inset_cells[[i]]$mar * canvas$canvas_scale
      } else {
        c(5.0, 5.0, 2.6, 1) * canvas$canvas_scale
      }
      old_mar <- par(mar = inset_mar)
    }
    panel <- if (single_locator_flags[i]) {
      render_network_glyph_panel(b$badge_rows, box_cols[i], value_label, not_value_label,
                                  network_scale, network_node_gap)
    } else {
      render_zoomed_scatter_panel(b$d_eligible, b$badge_rows, b$tie_box_rows, b$zbox, ranges, nofreq,
                                   axis_label, point_label, box_cols[i], badge_selection, badge_layout,
                                   value_label, not_value_label, neutral, max_group_members,
                                   zoom_annotated_only = zoom_annotated_only,
                                   repel_annotated_only = repel_annotated_only,
                                   kld_variant = kld_variant, label_lang = label_lang)
    }
    if (explicit && single_locator_flags[i] && !lay$inset_cells[[i]]$hide_x) {
      kld_sub <- if (nofreq) "nofreq" else "signed"
      title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"),
            line = par("mgp")[1] + 1.4)
      options(kld.plot.hide_y_title = TRUE)
      .draw_axis_direction_arrows(value_label, not_value_label)
    }
    if (multi) {
      .draw_inset_header(b$label, b$desc, box_cols[i], cex = inset_header_cex,
                         label_cex = inset_label_cex)
      par(old_mar)
    }
    if (length(panel$pending) > 0) {
      all_pending[[length(all_pending) + 1]] <- panel$pending
      all_marks_in[[length(all_marks_in) + 1]] <-
        .draw_svg_calibration_marks(panel$in_to_data_x, panel$in_to_data_y, panel$pin)
    }
    if (!multi) {
      single_dest_ndc_x <- panel$dest_ndc_x
      single_dest_ndc_y <- panel$dest_ndc_y
    }
  }

  if (!multi) {
    .draw_projection_lines(source_panel$src_ndc_x[[1]], source_panel$src_ndc_y[[1]],
                            single_dest_ndc_x, single_dest_ndc_y, inset_border_col)
  }

  dev.off()
  splice_ok <- TRUE
  if (length(all_pending) > 0) {
    .splice_network_glyphs(output, all_pending, all_marks_in)
  }
}

.kv_pairs <- function(tokens) {
  kv <- strsplit(tokens, "=", fixed = TRUE)
  setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))
}

.args_from_cli <- function() {
  .kv_pairs(commandArgs(trailingOnly = TRUE))
}

.inset_configs_from_args <- function(args) {
  keys <- names(args)
  m <- regmatches(keys, regexec("^inset([0-9]+)_(.+)$", keys))
  has_match <- vapply(m, function(x) length(x) == 3, logical(1))
  if (!any(has_match)) return(NULL)
  idx_num <- as.integer(vapply(m[has_match], `[`, "", 2))
  sub_keys <- vapply(m[has_match], `[`, "", 3)
  vals <- unname(args[has_match])
  by_idx <- split(seq_along(idx_num), idx_num)
  group_order <- names(by_idx)[order(as.integer(names(by_idx)))]
  lapply(group_order, function(gi) {
    rows <- by_idx[[gi]]
    get1 <- function(name) {
      hit <- which(sub_keys[rows] == name)
      if (length(hit) == 0) NULL else vals[rows][hit[1]]
    }
    filt_val <- function(type) {
      v <- get1(paste0(type, "_threshold")); if (is.null(v)) NULL else as.numeric(v)
    }
    filt_op <- function(type) {
      v <- get1(paste0(type, "_threshold_op")); if (is.null(v)) "gte" else v
    }
    list(
      inset_position = get1("inset_position"),
      frontier_n = { v <- get1("frontier_n"); if (is.null(v)) NULL else as.numeric(v) },
      filters = .selection_filter_spec(
        edge = filt_val("edge"), edge_op = filt_op("edge"),
        dispersion = filt_val("dispersion"), dispersion_op = filt_op("dispersion"),
        node_dispersion = filt_val("node_dispersion"), node_dispersion_op = filt_op("node_dispersion"),
        frequency = filt_val("frequency"), frequency_op = filt_op("frequency"),
        node_frequency = filt_val("node_frequency"), node_frequency_op = filt_op("node_frequency"),
        node_keyness = filt_val("node_keyness"), node_keyness_op = filt_op("node_keyness"),
        node_signed_keyness = filt_val("node_signed_keyness"),
        node_signed_keyness_op = filt_op("node_signed_keyness")
      )
    )
  })
}

.draw_point_vs_target_args <- function(args, output = args[["output"]]) {
  .arg <- function(name, default = NULL) {
    if (name %in% names(args)) args[[name]] else default
  }
  inset_configs <- .inset_configs_from_args(args)
  list(
    input = args[["input"]], output = output,
    point_label = .arg("point_label", "point"),
    axis_label = .arg("axis_label", "target"),
    value_label = .arg("value_label", "value"),
    not_value_label = .arg("not_value_label", paste0("¬", .arg("value_label", "value"))),
    nofreq = tolower(.arg("nofreq", "false")) == "true",
    badge_layout = .arg("badge_layout", "empty_region"),
    annotate_words_file = .arg("annotate_words_file"),
    highlight_file = .arg("highlight_file"),
    highlight_pairs_file = .arg("highlight_pairs_file"),
    annotate_pairs_file = .arg("annotate_pairs_file"),
    annotate_nontarget_only = tolower(.arg("annotate_nontarget_only", "false")) == "true",
    annotate_point_exclude = if (!is.null(.arg("annotate_point_exclude"))) {
      strsplit(.arg("annotate_point_exclude"), ",", fixed = TRUE)[[1]]
    } else NULL,
    highlight_frontier = tolower(.arg("highlight_frontier", "false")) == "true",
    default_label = .arg("default_label", .arg("point_label", "point")),
    badge_selection = .arg("badge_selection", "words"),
    frontier_n = as.numeric(.arg("frontier_n", "15")),
    inset = tolower(.arg("inset", "false")) == "true" || !is.null(inset_configs),
    inset_position = .arg("inset_position", "top"),
    inset_configs = inset_configs,
    inset_pad_frac = as.numeric(.arg("inset_pad_frac", "0.15")),
    zoom_annotated_only = tolower(.arg("zoom_annotated_only", "false")) == "true",
    repel_annotated_only = tolower(.arg("repel_annotated_only", "false")) == "true",
    edge_threshold = .arg("edge_threshold"),
    edge_threshold_op = .arg("edge_threshold_op", "gte"),
    dispersion_threshold = .arg("dispersion_threshold"),
    dispersion_threshold_op = .arg("dispersion_threshold_op", "gte"),
    node_dispersion_threshold = .arg("node_dispersion_threshold"),
    node_dispersion_threshold_op = .arg("node_dispersion_threshold_op", "gte"),
    frequency_threshold = .arg("frequency_threshold"),
    frequency_threshold_op = .arg("frequency_threshold_op", "gte"),
    node_frequency_threshold = .arg("node_frequency_threshold"),
    node_frequency_threshold_op = .arg("node_frequency_threshold_op", "gte"),
    node_keyness_threshold = .arg("node_keyness_threshold"),
    node_keyness_threshold_op = .arg("node_keyness_threshold_op", "gte"),
    node_signed_keyness_threshold = .arg("node_signed_keyness_threshold"),
    node_signed_keyness_threshold_op = .arg("node_signed_keyness_threshold_op", "gte"),
    aggregate_locators = tolower(.arg("aggregate_locators", "true")) == "true",
    max_group_members = as.numeric(.arg("max_group_members", "50")),
    network_scale = as.numeric(.arg("network_scale", "0.01")),
    network_node_gap = as.numeric(.arg("network_node_gap", "1")),
    point_field = .arg("point_field", "entry"),
    kld_variant = .arg("kld_variant"),
    frequency_type = .arg("frequency_type"),
    label_lang = .arg("label_lang", "ja"),
    unselected_alpha = as.numeric(.arg("unselected_alpha", "0.05")),
    main_position = .arg("main_position"),
    inset_header_cex = as.numeric(.arg("inset_header_cex", "0.75")),
    inset_label_cex = { v <- .arg("inset_label_cex"); if (is.null(v)) NULL else as.numeric(v) },
    box_label_cex = as.numeric(.arg("box_label_cex", "0.8")),
    text_scale = as.numeric(.arg("text_scale", "1")),
    hide_x_title = tolower(.arg("hide_x_title", "false")) == "true",
    hide_y_title = tolower(.arg("hide_y_title", "false")) == "true",
    hide_x_ticks = tolower(.arg("hide_x_ticks", "false")) == "true",
    hide_y_ticks = tolower(.arg("hide_y_ticks", "false")) == "true"
  )
}

if (sys.nframe() == 0 && !interactive()) {
  do.call(draw_point_vs_target, .draw_point_vs_target_args(.args_from_cli()))
}
