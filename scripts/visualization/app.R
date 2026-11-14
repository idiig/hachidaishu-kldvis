
library(shiny)

if (exists(".base_url", envir = .GlobalEnv, inherits = FALSE)) {
  downloadButton <- function(...) {
    tag <- shiny::downloadButton(...)
    tag$attribs$download <- NULL
    tag
  }
}

APP_DIR <- getwd()
ROOT_DIR <- normalizePath(file.path(APP_DIR, "..", ".."))

source(file.path(APP_DIR, "plot_point_vs_target.R"))
options(kld.plot.font_family = Sys.getenv("KLD_PLOT_FONT", "Noto Serif CJK JP"))
if (FALSE) {
  library(ggrepel)
  library(data.table)
  library(igraph)
  library(jsonlite)
}

VIS_DIR <- file.path(ROOT_DIR, "data", "processed", "visualization")
REF_DIR <- file.path(ROOT_DIR, "data", "reference")
.MANIFEST_PATH <- file.path(VIS_DIR, "manifest.json")
.DATA_MANIFEST <- if (file.exists(.MANIFEST_PATH)) jsonlite::read_json(.MANIFEST_PATH) else list()
.table_header <- function(path) {
  meta <- .DATA_MANIFEST[[basename(path)]]
  if (!is.null(meta)) return(unlist(meta$columns, use.names = FALSE))
  strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
}

.BLANK_PREVIEW_PNG <- local({
  path <- tempfile(fileext = ".png")
  png(path, width = 1, height = 1, bg = "transparent")
  dev.off()
  path
})

.NONE <- "(none)"

.annotate_words_choices <- function(dir) {
  c(.NONE, list.files(dir, pattern = "\\.txt$", full.names = FALSE))
}

.tsv_choices_by_header <- function(dir, required_cols) {
  files <- list.files(dir, pattern = "\\.tsv$", full.names = FALSE)
  matching <- Filter(function(f) {
    header <- strsplit(readLines(file.path(dir, f), n = 1, warn = FALSE), "\t")[[1]]
    all(required_cols %in% header)
  }, files)
  c(.NONE, matching)
}

.resolve_ref_choice <- function(dir, choice) {
  if (identical(choice, .NONE) || is.null(choice) || choice == "") return(NULL)
  file.path(dir, choice)
}

.kv_from_cmdline <- function(text) {
  tokens <- strsplit(trimws(text), "\\s+")[[1]]
  tokens <- tokens[nzchar(tokens) & grepl("=", tokens, fixed = TRUE)]
  if (length(tokens) == 0) return(setNames(character(0), character(0)))
  .kv_pairs(tokens)
}

.resolve_cmdline_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (startsWith(path, "/")) path else file.path(ROOT_DIR, path)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.sync_ui_from_kv <- function(session, kv) {
  .get <- function(key) if (key %in% names(kv)) kv[[key]] else NULL
  .bool <- function(key, default = FALSE) {
    v <- .get(key)
    if (is.null(v)) default else tolower(v) == "true"
  }
  .num <- function(key) {
    v <- .get(key)
    if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
  }
  .select_file <- function(id, key) {
    v <- .get(key)
    updateSelectInput(session, id, selected = if (is.null(v) || !nzchar(v)) .NONE else basename(v))
  }
  .sync_filter <- function(id) {
    v <- .num(id)
    if (!is.na(v)) {
      updateCheckboxInput(session, paste0("enable_", id), value = TRUE)
      updateNumericInput(session, id, value = v)
      op <- .get(paste0(id, "_op"))
      updateSelectInput(session, paste0(id, "_op"), selected = if (is.null(op)) "gte" else op)
    } else {
      updateCheckboxInput(session, paste0("enable_", id), value = FALSE)
    }
  }

  if (!is.null(.get("input"))) {
    updateSelectInput(session, "input_file", selected = basename(.get("input")))
  }
  if (!is.null(.get("kld_variant"))) {
    updateSelectInput(session, "kld_variant", selected = .get("kld_variant"))
  }
  if (!is.null(.get("frequency_type"))) {
    updateSelectInput(session, "frequency_type", selected = .get("frequency_type"))
  }
  updateTextInput(session, "point_label", value = .get("point_label") %||% "point")
  updateTextInput(session, "axis_label", value = .get("axis_label") %||% "target")
  updateTextInput(session, "value_label", value = .get("value_label") %||% "value")
  if (!is.null(.get("badge_selection"))) {
    updateSelectInput(session, "badge_selection", selected = .get("badge_selection"))
  }
  updateCheckboxInput(session, "aggregate_locators", value = .bool("aggregate_locators", TRUE))

  for (id in c("edge_threshold", "dispersion_threshold", "node_dispersion_threshold",
               "frequency_threshold", "node_frequency_threshold", "node_keyness_threshold")) {
    .sync_filter(id)
  }
  frontier_n <- .num("frontier_n")
  updateNumericInput(session, "frontier_n", value = if (is.na(frontier_n)) 15 else frontier_n)
  updateCheckboxInput(session, "annotate_nontarget_only", value = .bool("annotate_nontarget_only"))
  updateTextInput(session, "annotate_point_exclude", value = .get("annotate_point_exclude") %||% "")
  updateCheckboxInput(session, "highlight_frontier", value = .bool("highlight_frontier"))
  unselected_alpha <- .num("unselected_alpha")
  updateNumericInput(session, "unselected_alpha",
                      value = if (is.na(unselected_alpha)) 0.05 else unselected_alpha)

  .select_file("annotate_words_file", "annotate_words_file")
  .select_file("highlight_file", "highlight_file")
  .select_file("highlight_pairs_file", "highlight_pairs_file")
  .select_file("annotate_pairs_file", "annotate_pairs_file")

  inset_idx_strs <- unique(sub("^inset([0-9]+)_.*$", "\\1", grep("^inset[0-9]+_", names(kv), value = TRUE)))
  inset_idx <- as.integer(inset_idx_strs)
  updateCheckboxInput(session, "inset", value = .bool("inset") || length(inset_idx) > 0)
  if (!is.null(.get("inset_position"))) {
    updateSelectInput(session, "inset_position", selected = .get("inset_position"))
  }
  inset_pad_frac <- .num("inset_pad_frac")
  updateNumericInput(session, "inset_pad_frac", value = if (is.na(inset_pad_frac)) 0.15 else inset_pad_frac)
  updateCheckboxInput(session, "zoom_annotated_only", value = .bool("zoom_annotated_only"))
  updateCheckboxInput(session, "repel_annotated_only", value = .bool("repel_annotated_only"))

  if (length(inset_idx) > 0) {
    updateNumericInput(session, "inset_n", value = max(inset_idx))
    for (i in inset_idx) {
      .id <- function(suf) sprintf("inset%d_%s", i, suf)
      pos <- .get(.id("inset_position"))
      if (!is.null(pos)) updateSelectInput(session, .id("inset_position"), selected = pos)
      fn <- .num(.id("frontier_n"))
      if (!is.na(fn)) updateNumericInput(session, .id("frontier_n"), value = fn)
      for (type in c("edge", "dispersion", "node_dispersion", "frequency", "node_frequency",
                     "node_keyness")) {
        tid <- .id(paste0(type, "_threshold"))
        v <- .num(tid)
        if (!is.na(v)) {
          updateCheckboxInput(session, paste0("enable_", tid), value = TRUE)
          updateNumericInput(session, tid, value = v)
          op <- .get(paste0(tid, "_op"))
          updateSelectInput(session, paste0(tid, "_op"), selected = if (is.null(op)) "gte" else op)
        } else {
          updateCheckboxInput(session, paste0("enable_", tid), value = FALSE)
        }
      }
    }
  } else {
    updateNumericInput(session, "inset_n", value = 1)
  }
}

.DATASET_LABELS <- list(
  "entry-vs-anthology-shinkokin.tsv" = list(point_label = "entry", axis_label = "anthology", value_label = "shinkokin"),
  "entry-vs-gender-first4.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "entry-vs-gender-kokin.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "entry-vs-budatemajor-koi-first4.tsv" = list(point_label = "entry", axis_label = "budatemajor", value_label = "恋"),
  "budatemajor-vs-gender-first4.tsv" = list(point_label = "budatemajor", axis_label = "gender", value_label = "F"),
  "concept-mostnarrow-vs-gender-first4.tsv" = list(point_label = "concept", axis_label = "gender", value_label = "F"),
  "concept-mostnarrow-vs-anthology-shinkokin.tsv" = list(point_label = "concept", axis_label = "anthology", value_label = "shinkokin"),
  "pair-entry-entry-vs-anthology-shinkokin.tsv" = list(point_label = "entry", axis_label = "anthology", value_label = "shinkokin"),
  "pair-entry-entry-vs-gender-first4.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "pair-entry-budatemajor-vs-gender-first4.tsv" = list(point_label = "budatemajor", axis_label = "gender", value_label = "F")
)

.infer_labels <- function(input_file, root_dir = ROOT_DIR) {
  .DATASET_LABELS[[input_file]]
}

.DEFAULT_INPUT_FILE <- "entry-vs-anthology-shinkokin.tsv"
.INITIAL_LABELS <- .infer_labels(.DEFAULT_INPUT_FILE)
if (is.null(.INITIAL_LABELS)) {
  .INITIAL_LABELS <- list(point_label = "point", axis_label = "target", value_label = "value")
}

.plottable_tables <- function(dir) {
  files <- list.files(dir, pattern = "\\.tsv$", full.names = FALSE)
  files <- union(names(.DATA_MANIFEST), files)
  Filter(function(f) "point" %in% .table_header(file.path(dir, f)), files)
}

.is_unenriched_pair <- function(path) {
  header <- .table_header(path)
  "w1" %in% header && !("freq1" %in% header)
}

.pair_header_fields <- function(path) {
  .table_header(path)
}

.is_entry_entry_pair <- function(path, header) {
  if (!all(c("field1", "field2") %in% header)) return(TRUE)
  meta <- .DATA_MANIFEST[[basename(path)]]
  if (!is.null(meta)) return(identical(meta$field1, "entry") && identical(meta$field2, "entry"))
  row <- read.table(path, sep = "\t", header = TRUE, nrows = 1,
                     stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  identical(row$field1[1], "entry") && identical(row$field2[1], "entry")
}

.cmd_toggle <- function(key) {
  tags$div(style = "margin: -6px 0 10px 18px; font-size: 0.85em; color: #666;",
           checkboxInput(paste0("cmd_incl_", key), "再現コマンドに含める", value = TRUE))
}

.THRESHOLD_OP_CHOICES <- c("以上 (≥)" = "gte", "以下 (≤)" = "lte")

.filter_control <- function(id, label, pair_only = FALSE) {
  body <- tagList(
    checkboxInput(paste0("enable_", id), paste0(label, "を有効にする"), value = FALSE),
    conditionalPanel(
      condition = paste0("input.enable_", id),
      numericInput(id, label, value = NA),
      selectInput(paste0(id, "_op"), NULL, choices = .THRESHOLD_OP_CHOICES, selected = "gte"),
      .cmd_toggle(id),
      .cmd_toggle(paste0(id, "_op"))
    )
  )
  if (pair_only) conditionalPanel(condition = "output.is_pair_input", body) else body
}

.inset_config_filter_control <- function(id, label) {
  tagList(
    checkboxInput(paste0("enable_", id), paste0(label, "を有効にする"), value = FALSE),
    conditionalPanel(
      condition = paste0("input.enable_", id),
      fluidRow(
        column(7, numericInput(id, label, value = NA)),
        column(5, selectInput(paste0(id, "_op"), NULL, choices = .THRESHOLD_OP_CHOICES, selected = "gte"))
      )
    )
  )
}

.resolve_threshold <- function(input, id) {
  if (!isTRUE(input[[paste0("enable_", id)]])) return(NULL)
  value <- input[[id]]
  if (is.na(value)) NULL else value
}

.resolve_threshold_op_input <- function(input, id) {
  if (!isTRUE(input[[paste0("enable_", id)]])) return(NULL)
  v <- input[[paste0(id, "_op")]]
  if (is.null(v)) "gte" else v
}

.INSET_POSITIONS <- c("top", "top_right", "right", "bottom_right",
                      "bottom", "bottom_left", "left", "top_left")
.MAX_INSETS <- length(.INSET_POSITIONS)

.THRESHOLD_TYPES <- c("edge", "dispersion", "node_dispersion", "frequency", "node_frequency",
                      "node_keyness")
.threshold_key_pairs <- function() {
  as.vector(rbind(paste0(.THRESHOLD_TYPES, "_threshold"), paste0(.THRESHOLD_TYPES, "_threshold_op")))
}

.INSET_CONFIG_SUBKEYS <- c("inset_position", "frontier_n", .threshold_key_pairs())

.inset_config_keys <- function() {
  unlist(lapply(seq_len(.MAX_INSETS), function(i) sprintf("inset%d_%s", i, .INSET_CONFIG_SUBKEYS)))
}

ALL_PARAM_KEYS <- c(
  "kld_variant", "frequency_type", "point_label", "axis_label", "value_label", "badge_selection",
  "aggregate_locators", .threshold_key_pairs(), "frontier_n", "annotate_nontarget_only",
  "annotate_point_exclude", "highlight_frontier", "unselected_alpha",
  "annotate_words_file", "highlight_file",
  "highlight_pairs_file", "annotate_pairs_file",
  "inset", "inset_position", "inset_pad_frac",
  "zoom_annotated_only", "repel_annotated_only",
  .inset_config_keys()
)

.build_repro_cmd <- function(input_path, params, excluded_keys) {
  kv <- Filter(Negate(is.null), lapply(ALL_PARAM_KEYS, function(key) {
    if (key %in% excluded_keys) return(NULL)
    value <- params[[key]]
    if (is.null(value) || identical(value, "") || identical(value, FALSE)) return(NULL)
    sprintf("%s=%s", key, paste(value, collapse = ","))
  }))
  paste(c(
    "Rscript scripts/visualization/plot_point_vs_target.R",
    sprintf("input=%s", input_path),
    unlist(kv)
  ), collapse = " ")
}

.repro_cmd_visual_width <- function(s) {
  sum(ifelse(utf8ToInt(enc2utf8(s)) > 0x7F, 2, 1))
}

.append_repro_cmd_line <- function(svg_text, cmd, row_pt = 22) {
  if (!nzchar(cmd)) return(svg_text)
  pattern <- paste0(
    '(<svg[^>]*width="([0-9.]+)pt"[^>]*height=")([0-9.]+)',
    '(pt"[^>]*viewBox="[0-9.\\-]+ [0-9.\\-]+ [0-9.\\-]+ )([0-9.]+)(")'
  )
  m <- regmatches(svg_text, regexec(pattern, svg_text))[[1]]
  if (length(m) == 0) return(svg_text)
  width <- as.numeric(m[3])
  old_height <- as.numeric(m[4])
  new_height <- old_height + row_pt
  old_root <- m[1]
  new_root <- paste0(m[2], new_height, m[5], new_height, m[7])
  svg_text <- sub(old_root, new_root, svg_text, fixed = TRUE)
  font_size <- width / (max(1, .repro_cmd_visual_width(cmd)) * 0.6)
  text_el <- sprintf(
    '<text x="2" y="%f" font-family="monospace" font-size="%f">%s</text></svg>',
    old_height + row_pt * 0.7, font_size, htmltools::htmlEscape(cmd)
  )
  sub("</svg>\\s*$", text_el, svg_text)
}

ui <- fluidPage(
  tags$head(tags$style(HTML(
    "#generating-msg { display: none; color: #c00; font-weight: bold; }
     html.shiny-busy #generating-msg { display: block; }"
  ))),
  titlePanel("point-vs-target 調整パネル"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("input_mode", "入力方式",
                   choices = c("UI" = "ui", "コマンドライン" = "cmdline"), selected = "ui",
                   inline = TRUE),
      conditionalPanel(
        condition = "input.input_mode == 'ui'",
      selectInput("input_file", "入力ファイル (data/processed/visualization/)",
                  choices = .plottable_tables(VIS_DIR),
                  selected = .DEFAULT_INPUT_FILE),
      uiOutput("unenriched_warning"),
      selectInput("kld_variant", "kld_variant", choices = c("obs", "nofreq", "norm", "keyassoc"),
                  selected = "keyassoc"),
      .cmd_toggle("kld_variant"),
      selectInput("frequency_type", "frequency_type", choices = c("term", "doc"),
                  selected = "term"),
      .cmd_toggle("frequency_type"),
      textInput("point_label", "point_label", .INITIAL_LABELS$point_label),
      .cmd_toggle("point_label"),
      textInput("axis_label", "axis_label", .INITIAL_LABELS$axis_label),
      .cmd_toggle("axis_label"),
      textInput("value_label", "value_label", .INITIAL_LABELS$value_label),
      .cmd_toggle("value_label"),
      selectInput("badge_selection", "badge_selection",
                  choices = c("words", "frontier", "pairs"), selected = "words"),
      .cmd_toggle("badge_selection"),
      checkboxInput("aggregate_locators", "aggregate_locators", value = TRUE),
      .cmd_toggle("aggregate_locators"),
      tags$hr(),
      tags$b("フィルタ (filter)"),
      .filter_control("edge_threshold", "edge_threshold (pair only)", pair_only = TRUE),
      .filter_control("dispersion_threshold", "dispersion_threshold"),
      .filter_control("node_dispersion_threshold", "node_dispersion_threshold (pair片側)",
                      pair_only = TRUE),
      .filter_control("frequency_threshold", "frequency_threshold (total_freq)"),
      .filter_control("node_frequency_threshold", "node_frequency_threshold (pair片側freq1/freq2)",
                      pair_only = TRUE),
      .filter_control("node_keyness_threshold",
                      "node_keyness_threshold (pair片側語自身のKLD、絶対値)",
                      pair_only = TRUE),
      conditionalPanel(
        condition = "input.badge_selection == 'frontier'",
        numericInput("frontier_n", "frontier_n (badge_selection=\"frontier\"時)", value = 15,
                     min = 1, step = 1),
        .cmd_toggle("frontier_n")
      ),
      checkboxInput("highlight_frontier",
                    "highlight_frontier (frequencyがfrontier中点より大きい点を染色、badge_selectionと独立)",
                    value = FALSE),
      .cmd_toggle("highlight_frontier"),
      numericInput("unselected_alpha",
                   "unselected_alpha (highlight/highlight_frontier選択時、非選択点の透明度)",
                   value = 0.05, min = 0, max = 1, step = 0.05),
      .cmd_toggle("unselected_alpha"),
      conditionalPanel(
        condition = "input.badge_selection == 'words' && input.annotate_words_file != '(none)'",
        checkboxInput("annotate_nontarget_only",
                      "annotate_nontarget_only (n_target==0の候補のみbadge化)", value = FALSE),
        .cmd_toggle("annotate_nontarget_only"),
        textInput("annotate_point_exclude",
                 "annotate_point_exclude (除外するpoint文字列、カンマ区切り)", value = ""),
        .cmd_toggle("annotate_point_exclude")
      ),
      tags$hr(),
      tags$b("アノテーションリスト (data/reference/)"),
      conditionalPanel(
        condition = "input.badge_selection == 'words'",
        selectInput("annotate_words_file", "annotate_words_file",
                    choices = .annotate_words_choices(REF_DIR)),
        .cmd_toggle("annotate_words_file")
      ),
      conditionalPanel(
        condition = "!output.is_pair_input",
        selectInput("highlight_file", "highlight_file",
                    choices = .tsv_choices_by_header(REF_DIR, c("term", "type"))),
        .cmd_toggle("highlight_file")
      ),
      conditionalPanel(
        condition = "output.is_entry_entry_pair_input",
        selectInput("highlight_pairs_file", "highlight_pairs_file",
                    choices = .tsv_choices_by_header(REF_DIR, c("w1", "w2", "type"))),
        .cmd_toggle("highlight_pairs_file")
      ),
      conditionalPanel(
        condition = "input.badge_selection == 'pairs'",
        selectInput("annotate_pairs_file", "annotate_pairs_file (badge_selection=\"pairs\"時)",
                    choices = .tsv_choices_by_header(REF_DIR, c("w1", "w2"))),
        .cmd_toggle("annotate_pairs_file")
      ),
      tags$hr(),
      checkboxInput("inset", "拡大表示 (inset)", value = FALSE),
      .cmd_toggle("inset"),
      conditionalPanel(
        condition = "input.inset",
        conditionalPanel(
          condition = "input.inset_n < 2",
          selectInput("inset_position", "inset_position",
                      choices = .INSET_POSITIONS, selected = "top"),
          .cmd_toggle("inset_position")
        ),
        numericInput("inset_pad_frac", "inset_pad_frac", value = 0.15, step = 0.05),
        .cmd_toggle("inset_pad_frac"),
        checkboxInput("zoom_annotated_only", "zoom_annotated_only (拡大側は注釈された点のみ)",
                      value = FALSE),
        .cmd_toggle("zoom_annotated_only"),
        checkboxInput("repel_annotated_only", "repel_annotated_only (repel対象も注釈された点のみ)",
                      value = FALSE),
        .cmd_toggle("repel_annotated_only"),
        tags$hr(),
        tags$b("複数inset設定 (inset_configs)"),
        tags$p(style = "font-size: 0.85em; color: #666;",
               "inset個数を2以上にすると、各insetごとに別々のinset_position/",
               "frontier_n/filterを設定できます(主図は常に完全な点群を表示、",
               "絞り込みは各insetの候補選定にのみ適用されます——edge_threshold等の",
               "上のフィルタ欄は個数1の時だけ使われます)。"),
        numericInput("inset_n", "inset個数", value = 1, min = 1, max = .MAX_INSETS, step = 1),
        .cmd_toggle("inset_configs"),
        conditionalPanel(
          condition = "input.inset_n >= 2",
          lapply(seq_len(.MAX_INSETS), function(i) {
            conditionalPanel(
              condition = sprintf("input.inset_n >= %d", i),
              wellPanel(
                tags$b(sprintf("inset %d (%s)", i, LETTERS[i])),
                selectInput(sprintf("inset%d_inset_position", i), "inset_position",
                            choices = .INSET_POSITIONS,
                            selected = .INSET_POSITIONS[((i - 1) %% length(.INSET_POSITIONS)) + 1]),
                numericInput(sprintf("inset%d_frontier_n", i), "frontier_n", value = 15,
                             min = 1, step = 1),
                .inset_config_filter_control(sprintf("inset%d_edge_threshold", i), "edge_threshold"),
                .inset_config_filter_control(sprintf("inset%d_dispersion_threshold", i), "dispersion_threshold"),
                .inset_config_filter_control(sprintf("inset%d_node_dispersion_threshold", i), "node_dispersion_threshold"),
                .inset_config_filter_control(sprintf("inset%d_frequency_threshold", i), "frequency_threshold"),
                .inset_config_filter_control(sprintf("inset%d_node_frequency_threshold", i), "node_frequency_threshold"),
                .inset_config_filter_control(sprintf("inset%d_node_keyness_threshold", i), "node_keyness_threshold")
              )
            )
          })
        )
      )
      ),
      conditionalPanel(
        condition = "input.input_mode == 'cmdline'",
        tags$p(style = "font-size: 0.85em; color: #666;",
               "「再現コマンドをコピー」で出力されるのと同じ形式のkey=value",
               "を空白区切りで貼り付けてください(先頭の\"Rscript ...\"部分や",
               "output=はあっても無視されます、input=はdata/processed/",
               "visualization/配下のファイルをリポジトリルート相対で指定)。"),
        textAreaInput("cmdline_text", NULL, rows = 8, width = "100%",
                      placeholder = paste(
                        "input=data/processed/visualization/entry-vs-anthology-shinkokin.tsv",
                        "kld_variant=keyassoc frequency_type=term badge_selection=words"
                      ))
      ),
      tags$hr(),
      actionButton("render", "描画", class = "btn-primary"),
      br(), br(),
      downloadButton("download_svg", "SVGを保存"),
      br(), br(),
      actionButton("copy_cmd", "再現コマンドをコピー")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("プロット",
          tags$div(id = "generating-msg", "画像生成中..."),
          imageOutput("preview", height = "auto"),
          verbatimTextOutput("status")
        ),
        tabPanel("表", tableOutput("annotated_table"))
      )
    )
  ),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('copy_repro_cmd', function(cmd) {
      navigator.clipboard.writeText(cmd);
    });
  "))
)

server <- function(input, output, session) {
  output$unenriched_warning <- renderUI({
    path <- file.path(VIS_DIR, input$input_file)
    if (.is_unenriched_pair(path)) {
      tags$p(style = "color: darkorange;",
             "このデータにはラベル表示に必要な列が不足しています。",
             "エラーにはなりませんが、badge(ラベル)が1つも出ません",
             "必要な列を含むデータを選んでください。")
    }
  })

  output$is_pair_input <- reactive({
    "w1" %in% .pair_header_fields(file.path(VIS_DIR, input$input_file))
  })
  outputOptions(output, "is_pair_input", suspendWhenHidden = FALSE)

  output$is_entry_entry_pair_input <- reactive({
    path <- file.path(VIS_DIR, input$input_file)
    header <- .pair_header_fields(path)
    "w1" %in% header && .is_entry_entry_pair(path, header)
  })
  outputOptions(output, "is_entry_entry_pair_input", suspendWhenHidden = FALSE)

  observeEvent(input$input_file, {
    labels <- .infer_labels(input$input_file)
    if (is.null(labels)) return(invisible(NULL))
    updateTextInput(session, "point_label", value = labels$point_label)
    updateTextInput(session, "axis_label", value = labels$axis_label)
    updateTextInput(session, "value_label", value = labels$value_label)
  }, ignoreInit = TRUE)

  last_draw_args <- reactiveVal(NULL)
  status <- reactiveVal("")
  repro_cmd <- reactiveVal("")
  annotated_table <- reactiveVal(NULL)

  observeEvent(input$render, {
    if (!identical(input$input_mode, "cmdline")) return(invisible(NULL))
    .sync_ui_from_kv(session, .kv_from_cmdline(input$cmdline_text))
  })

  rendered <- eventReactive(input$render, {
    out <- tempfile(fileext = ".png")

    if (identical(input$input_mode, "cmdline")) {
      kv <- .kv_from_cmdline(input$cmdline_text)
      if (length(kv) == 0 || !("input" %in% names(kv))) {
        status("エラー: input=... を含むコマンドを貼り付けてください")
        return(NULL)
      }
      kv[["input"]] <- .resolve_cmdline_path(kv[["input"]])
      for (path_key in c("annotate_words_file", "highlight_file",
                          "highlight_pairs_file", "annotate_pairs_file")) {
        if (path_key %in% names(kv)) kv[[path_key]] <- .resolve_cmdline_path(kv[[path_key]])
      }
      draw_args <- .draw_point_vs_target_args(kv, output = out)
      ok <- tryCatch({
        annotated_table(do.call(draw_point_vs_target, draw_args))
        TRUE
      }, error = function(e) {
        status(paste("エラー:", conditionMessage(e)))
        FALSE
      })
      if (!ok || !file.exists(out)) return(NULL)
      status(paste("描画完了 (cmdline):", basename(kv[["input"]])))
      last_draw_args(draw_args)
      repro_cmd(paste(
        "Rscript scripts/visualization/plot_point_vs_target.R",
        paste(names(kv), kv, sep = "=", collapse = " ")
      ))
      return(out)
    }

    input_path <- file.path(VIS_DIR, input$input_file)
    edge_threshold <- .resolve_threshold(input, "edge_threshold")
    edge_threshold_op <- .resolve_threshold_op_input(input, "edge_threshold")
    dispersion_threshold <- .resolve_threshold(input, "dispersion_threshold")
    dispersion_threshold_op <- .resolve_threshold_op_input(input, "dispersion_threshold")
    node_dispersion_threshold <- .resolve_threshold(input, "node_dispersion_threshold")
    node_dispersion_threshold_op <- .resolve_threshold_op_input(input, "node_dispersion_threshold")
    frequency_threshold <- .resolve_threshold(input, "frequency_threshold")
    frequency_threshold_op <- .resolve_threshold_op_input(input, "frequency_threshold")
    node_frequency_threshold <- .resolve_threshold(input, "node_frequency_threshold")
    node_frequency_threshold_op <- .resolve_threshold_op_input(input, "node_frequency_threshold")
    node_keyness_threshold <- .resolve_threshold(input, "node_keyness_threshold")
    node_keyness_threshold_op <- .resolve_threshold_op_input(input, "node_keyness_threshold")
    inset_n <- input$inset_n
    inset_n <- if (is.null(inset_n) || is.na(inset_n)) 1 else inset_n
    inset_configs <- if (isTRUE(input$inset) && inset_n >= 2) {
      lapply(seq_len(inset_n), function(i) {
        .id <- function(suf) sprintf("inset%d_%s", i, suf)
        list(
          inset_position = input[[.id("inset_position")]],
          frontier_n = input[[.id("frontier_n")]],
          filters = setNames(lapply(c("edge", "dispersion", "node_dispersion",
                                       "frequency", "node_frequency", "node_keyness"), function(type) {
            tid <- .id(paste0(type, "_threshold"))
            list(value = .resolve_threshold(input, tid), op = .resolve_threshold_op_input(input, tid))
          }), c("edge", "dispersion", "node_dispersion", "frequency", "node_frequency", "node_keyness"))
        )
      })
    } else NULL
    exclude_raw <- trimws(input$annotate_point_exclude)
    annotate_point_exclude <- if (exclude_raw == "") NULL else strsplit(exclude_raw, ",", fixed = TRUE)[[1]]
    annotate_words_file <- .resolve_ref_choice(REF_DIR, input$annotate_words_file)
    highlight_file <- .resolve_ref_choice(REF_DIR, input$highlight_file)
    highlight_pairs_file <- .resolve_ref_choice(REF_DIR, input$highlight_pairs_file)
    annotate_pairs_file <- .resolve_ref_choice(REF_DIR, input$annotate_pairs_file)
    draw_args <- list(
      input = input_path, output = out,
      kld_variant = input$kld_variant,
      frequency_type = input$frequency_type,
      point_label = input$point_label, axis_label = input$axis_label,
      value_label = input$value_label,
      badge_selection = input$badge_selection,
      aggregate_locators = input$aggregate_locators,
      edge_threshold = edge_threshold, edge_threshold_op = edge_threshold_op,
      dispersion_threshold = dispersion_threshold, dispersion_threshold_op = dispersion_threshold_op,
      node_dispersion_threshold = node_dispersion_threshold,
      node_dispersion_threshold_op = node_dispersion_threshold_op,
      frequency_threshold = frequency_threshold, frequency_threshold_op = frequency_threshold_op,
      node_frequency_threshold = node_frequency_threshold,
      node_frequency_threshold_op = node_frequency_threshold_op,
      node_keyness_threshold = node_keyness_threshold,
      node_keyness_threshold_op = node_keyness_threshold_op,
      frontier_n = input$frontier_n,
      annotate_nontarget_only = input$annotate_nontarget_only,
      annotate_point_exclude = annotate_point_exclude,
      highlight_frontier = input$highlight_frontier,
      unselected_alpha = input$unselected_alpha,
      annotate_words_file = annotate_words_file,
      highlight_file = highlight_file,
      highlight_pairs_file = highlight_pairs_file,
      annotate_pairs_file = annotate_pairs_file,
      inset = input$inset, inset_position = input$inset_position,
      inset_configs = inset_configs,
      inset_pad_frac = input$inset_pad_frac,
      zoom_annotated_only = input$zoom_annotated_only,
      repel_annotated_only = input$repel_annotated_only
    )
    ok <- tryCatch({
      annotated_table(do.call(draw_point_vs_target, draw_args))
      TRUE
    }, error = function(e) {
      status(paste("エラー:", conditionMessage(e)))
      FALSE
    })
    if (!ok || !file.exists(out)) return(NULL)
    status(paste("描画完了:", input$input_file))
    last_draw_args(draw_args)

    .rel_ref <- function(choice) if (identical(choice, .NONE)) NULL else file.path("data", "reference", choice)
    cmd_params <- list(
      kld_variant = input$kld_variant,
      frequency_type = input$frequency_type,
      point_label = input$point_label, axis_label = input$axis_label,
      value_label = input$value_label, badge_selection = input$badge_selection,
      aggregate_locators = input$aggregate_locators,
      edge_threshold = edge_threshold, edge_threshold_op = edge_threshold_op,
      dispersion_threshold = dispersion_threshold, dispersion_threshold_op = dispersion_threshold_op,
      node_dispersion_threshold = node_dispersion_threshold,
      node_dispersion_threshold_op = node_dispersion_threshold_op,
      frequency_threshold = frequency_threshold, frequency_threshold_op = frequency_threshold_op,
      node_frequency_threshold = node_frequency_threshold,
      node_frequency_threshold_op = node_frequency_threshold_op,
      node_keyness_threshold = node_keyness_threshold,
      node_keyness_threshold_op = node_keyness_threshold_op,
      frontier_n = input$frontier_n,
      annotate_nontarget_only = input$annotate_nontarget_only,
      annotate_point_exclude = annotate_point_exclude,
      highlight_frontier = input$highlight_frontier,
      unselected_alpha = input$unselected_alpha,
      annotate_words_file = .rel_ref(input$annotate_words_file),
      highlight_file = .rel_ref(input$highlight_file),
      highlight_pairs_file = .rel_ref(input$highlight_pairs_file),
      annotate_pairs_file = .rel_ref(input$annotate_pairs_file),
      inset = input$inset, inset_position = input$inset_position,
      inset_pad_frac = input$inset_pad_frac,
      zoom_annotated_only = input$zoom_annotated_only,
      repel_annotated_only = input$repel_annotated_only
    )
    if (!is.null(inset_configs)) {
      for (i in seq_along(inset_configs)) {
        .id <- function(suf) sprintf("inset%d_%s", i, suf)
        cfg <- inset_configs[[i]]
        cmd_params[[.id("inset_position")]] <- cfg$inset_position
        cmd_params[[.id("frontier_n")]] <- cfg$frontier_n
        for (type in names(cfg$filters)) {
          tid <- .id(paste0(type, "_threshold"))
          cmd_params[[tid]] <- cfg$filters[[type]]$value
          cmd_params[[paste0(tid, "_op")]] <- cfg$filters[[type]]$op
        }
      }
    }
    excluded_keys <- Filter(
      function(k) {
        if (grepl("^inset[0-9]+_", k)) isFALSE(input[["cmd_incl_inset_configs"]])
        else isFALSE(input[[paste0("cmd_incl_", k)]])
      },
      ALL_PARAM_KEYS
    )
    cmd <- .build_repro_cmd(
      file.path("data", "processed", "visualization", input$input_file),
      cmd_params, excluded_keys
    )
    repro_cmd(cmd)
    out
  })

  output$preview <- renderImage({
    path <- rendered()
    if (is.null(path)) return(list(src = .BLANK_PREVIEW_PNG, contentType = "image/png"))
    list(src = path, contentType = "image/png", style = "max-width: 100%; height: auto;")
  }, deleteFile = FALSE)

  output$annotated_table <- renderTable({
    tbl <- annotated_table()
    if (is.null(tbl) || nrow(tbl) == 0) return(data.frame())
    tbl
  })

  output$status <- renderText(status())

  output$download_svg <- downloadHandler(
    filename = function() paste0(tools::file_path_sans_ext(input$input_file), ".svg"),
    content = function(file) {
      draw_args <- last_draw_args()
      if (is.null(draw_args)) stop("先に「描画」を押してください")
      svg_out <- tempfile(fileext = ".svg")
      draw_args$output <- svg_out
      do.call(draw_point_vs_target, draw_args)
      writeLines(
        .append_repro_cmd_line(paste(readLines(svg_out, warn = FALSE), collapse = "\n"), repro_cmd()),
        svg_out
      )
      file.copy(svg_out, file)
    }
  )

  observeEvent(input$copy_cmd, {
    session$sendCustomMessage("copy_repro_cmd", repro_cmd())
  })
}

shinyApp(ui, server)
