#!/usr/bin/env Rscript
#
# plot_theme.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Defines:
#   bc_palette     : named Okabe-Ito colourblind-safe palette
#   model_colours  : stable model -> colour map (clinical / EEG / combined)
#   theme_neo()    : theme_classic base + top legend + bold facet strips
#   number_ticks() : pretty()-based axis-tick helper
#   save_fig()     : writes BOTH a 600-dpi sRGB PNG and a vector PDF
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

# ---- Okabe-Ito palette (colourblind-safe). Preference order for series:
#   orange (1st) -> skyblue (2nd) -> green (3rd) -> blue (4th) -> vermillion ...
bc_palette <- c(grey = "#999999", orange = "#E69F00", skyblue = "#56B4E9",
                green = "#009E73", blue = "#0072B2", vermillion = "#D55E00",
                purple = "#CC79A7", yellow = "#F0E442")

# ---- Stable model colour mapping, reused by every figure that shows models.
model_colours <- c(
  "Clinical"        = unname(bc_palette["grey"]),
  "EEG"             = unname(bc_palette["orange"]),
  "Clinical + EEG"  = unname(bc_palette["blue"])
)

# ---- Shared theme: theme_classic + top legend + bold facet strips.
theme_neo <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      strip.background = ggplot2::element_blank(),
      strip.text.x    = ggplot2::element_text(face = "bold"),
      strip.text.y    = ggplot2::element_text(face = "bold"))
}

# ---- pretty()-based axis ticks (define once, reuse).
number_ticks <- function(n) function(limits) pretty(limits, n)

# ---- Probe whether cairo_pdf actually works at runtime. capabilities("cairo")
# can report TRUE while the cairo DLL fails to load (missing X11 libs on some
# macOS installs), so we test it for real once and cache the result.
.neo_cairo_ok <- local({
  cached <- NA
  function() {
    if (!is.na(cached)) return(cached)
    ok <- FALSE
    if (isTRUE(suppressWarnings(capabilities("cairo")))) {
      tf <- tempfile(fileext = ".pdf")
      ok <- tryCatch(
        withCallingHandlers({
          grDevices::cairo_pdf(tf)
          try(grDevices::dev.off(), silent = TRUE)
          isTRUE(file.exists(tf) && file.info(tf)$size > 0)
        }, warning = function(w) invokeRestart("muffleWarning")),
        error = function(e) FALSE)
      if (file.exists(tf)) unlink(tf)
    }
    cached <<- isTRUE(ok)
    if (!cached) message("[plot_theme] cairo_pdf unavailable; using base 'pdf' device + Ghostscript font embedding")
    cached
  }
})

# ---- Embed all fonts into a vector PDF, in place (base 'pdf' device does not).
.embed_pdf_fonts <- function(path) {
  cand <- c(Sys.getenv("R_GSCMD", ""), Sys.which("gs"), "/opt/homebrew/bin/gs", "/usr/local/bin/gs", "/usr/bin/gs")
  gs <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(gs)) { warning("[plot_theme] Ghostscript not found; fonts NOT embedded in ", basename(path), call. = FALSE); return(invisible(FALSE)) }
  tmp <- paste0(path, ".embed.tmp")
  status <- system2(gs, c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", "-dEmbedAllFonts=true",
                          "-dSubsetFonts=true", "-dPDFSETTINGS=/prepress", "-dCompatibilityLevel=1.4",
                          paste0("-sOutputFile=", shQuote(tmp)), shQuote(path)), stdout = FALSE, stderr = FALSE)
  if (status == 0 && file.exists(tmp) && file.info(tmp)$size > 0) file.rename(tmp, path) else if (file.exists(tmp)) unlink(tmp)
  invisible(status == 0)
}

# ---- Save a figure as BOTH a 600-dpi PNG (sRGB via ragg) and a vector PDF.
save_fig <- function(plot, name, width, height, output_dir = CFG$dir_fig) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(output_dir, paste0(name, ".png"))
  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot, device = ragg::agg_png, width = width,
                    height = height, dpi = 600, bg = "white")
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 600, bg = "white")
  }
  cairo <- .neo_cairo_ok()
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  device = if (cairo) grDevices::cairo_pdf else "pdf", bg = "white")
  if (!cairo) .embed_pdf_fonts(pdf_path)
  invisible(c(png_path, pdf_path))
}
