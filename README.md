# KLD visualization

Live plotting app: **https://idiig.github.io/kld-visualization/**

Source: https://github.com/idiig/kld-visualization

A standalone R plotting and online parameter-tuning app. GitHub Pages runs R in the browser via Shinylive/webR; local and Docker use plain Shiny. Both modes share the same plotting code.

This repository contains only the display code, precomputed plotting data, and annotation lists — it does not include the statistical calculation pipeline, raw corpus, indices, or calculation configs.

## Nix

```sh
nix develop
Rscript tests/smoke.R
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Open http://localhost:3838 . If flakes aren't enabled on your Nix install, append `--extra-experimental-features 'nix-command flakes'` after `nix`.

To build and preview the same browser build used on GitHub Pages:

```sh
nix develop --command Rscript scripts/build_site.R
nix develop --command Rscript -e 'httpuv::runStaticServer("_site", port=8080)'
```

Open http://localhost:8080 . The first build needs internet access to download Shinylive and WebAssembly packages.

## Docker

```sh
docker compose up --build app
```

Open http://localhost:3838 . To build and preview the browser version directly:

```sh
docker compose --profile pages up --build pages
```

Open http://localhost:8080 . Docker and Nix are local/build environments only; GitHub Pages itself does not run a container or an R server.

## GitHub Pages

1. Push this directory as a standalone Git repository to the `main` branch on GitHub.
2. In the repository, go to **Settings → Pages → Build and deployment → Source** and select **GitHub Actions**.
3. `.github/workflows/pages.yaml` automatically runs the plotting checks, exports `_site`, and publishes it. Pull requests only run the build check.

The live URL is given by the Actions deployment output. No fixed repository path needs to be set, and project-Pages subpaths are supported.

## Data and plotting

`data/processed/visualization/*.rds` are lossless serialized versions of the precomputed TSVs, keeping all rows and columns; `manifest.json` records row counts, column names, and the MD5 of the source TSV. The app keeps the original TSV names, so existing commands can still read the corresponding RDS:

```sh
Rscript scripts/visualization/plot_point_vs_target.R \
  input=data/processed/visualization/entry-vs-gender-kokin.tsv \
  output=artifacts/example.svg kld_variant=norm frequency_type=term
```

To just repackage externally precomputed TSVs (without running any statistical calculation):

```sh
Rscript scripts/package_data.R /path/to/precomputed/visualization
```

`data/reference/` holds the highlight and annotation lists. Automatic labels are display metadata kept independently inside the app.

The browser version keeps parameter adjustment, filtering, zoomed insets, network diagrams, result tables, reproduction commands, and SVG download. Network diagrams are drawn entirely with R/igraph; the layout differs from the original Graphviz version and no external `neato` is needed. Japanese PNG/SVG output uses the bundled Noto Serif CJK JP font (license in `fonts/OFL.txt`); the font for local environments can be configured via `KLD_PLOT_FONT`.

The browser downloads each selected table on demand, and downloaded tables are kept for the current browser session. Large pair-word tables use more memory, so a desktop browser is recommended. Plotting is triggered by the "描画" (render) button.

## Browser verification

The native plotting check is `Rscript tests/smoke.R`. End-to-end verification uses Chromium:

```sh
npm install --no-save playwright
npx playwright install chromium
KLD_SERVE=1 node tests/browser.mjs http://127.0.0.1:8080/
```

`KLD_SERVE=1` automatically starts/stops a local static server for `_site` and requires Python 3; alternatively, host `_site` yourself and pass the URL to the test script. The test covers startup, dataset selection, PNG plotting, re-render on parameter change, SVG download, the result table, and word-pair plotting.
