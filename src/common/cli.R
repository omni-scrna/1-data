# Shared CLI helpers for omnibenchmark module entrypoints (R).
#
# Reserved, overwrite-on-update path (src/common/r/) — see AGENTS.md. The author
# owns an `argparser` parser; these helpers just add the benchmark's shared flags
# (universal base + the stage's I/O contract, from JSON in src/common/schema/) onto
# it: each schema arg becomes one argparser::add_argument call, then argparser owns
# naming, typing and parsing (no required/choices/dest handling — that's the Python
# engine). The author adds their own params and parses directly:
#
#   cli <- new.env(); source("src/common/cli.R", local = cli)
#   p <- arg_parser("PCA module")
#   p <- cli$add_base_args(p)                  # --output_dir, --name
#   p <- cli$add_stage_args(p, "embedding")    # the stage I/O contract
#   p <- add_argument(p, "--n_components", type = "integer", help = "PCs")  # your own
#   args <- parse_args(p)                       # argparser's own parser
#
# These files are *copied* in by scripts/pull.py, not hand-written. See
# docs/common-code.md and docs/cli.md. Schema arg-spec: {flag, type, help?};
# types path|string|integer|number (path/string -> character).

suppressPackageStartupMessages(library(argparser))

COMMON_VERSION <- "0.1.0"  # x-release-version — stamped from src/common/VERSION by `pixi run version`
SCHEMA_DIR <- "src/common/schema"   # copied-in layout; override per call or set once after sourcing

common_version <- function() COMMON_VERSION

`%||%` <- function(a, b) if (is.null(a)) b else a

# Add a schema file's args onto the author's parser; argparser does the rest.
.inject <- function(p, path) {
  if (!file.exists(path)) stop(sprintf("schema not found: %s", path), call. = FALSE)
  for (a in jsonlite::fromJSON(path, simplifyDataFrame = FALSE)$args)
    p <- argparser::add_argument(p, a$flag, help = a$help %||% "",
           type = switch(a$type %||% "string", integer = "integer", number = "numeric", "character"))
  p
}

# On purpose, we don't try to enforce options or choices — left for later.

add_base_args  <- function(p, schema_dir = SCHEMA_DIR) .inject(p, file.path(schema_dir, "_base.json"))
add_stage_args <- function(p, interface, schema_dir = SCHEMA_DIR) .inject(p, file.path(schema_dir, paste0(interface, ".json")))
