# Shared CLI helpers for omnibenchmark module entrypoints (R).
#
# Reserved, overwrite-on-update path (src/common/r/) — see AGENTS.md. The author
# owns an `argparser` parser; these helpers just add the benchmark's shared flags
# (universal base + the stage's I/O contract, from JSON in src/common/schema/) onto
# it: each schema arg becomes one argparser::add_argument call, then argparser owns
# naming, typing and parsing (no required/choices/dest handling — that's the Python
# engine). The author adds their own params and parses directly:
#
# source("src/common/cli.R")
# p <- arg_parser("my module")
# p <- add_base_args(p)                    # --output_dir, --name
# p <- add_stage_args(p, "two-filter")     # the stage I/O contract
## your own method params — argparser directly (its add_argument requires `help`):
# p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
# args <- parse_args(p)
#
# These files are *copied* in by scripts/pull.py, not hand-written. See
# docs/common-code.md and docs/cli.md.

suppressPackageStartupMessages(library(argparser))

# where the schemas were copied on sync
SCHEMA_DIR <- "src/common/schema"

# Add a schema file's args onto the author's parser; argparser does the rest.
.inject <- function(p, path) {
  if (!file.exists(path)) stop(sprintf("schema not found: %s", path), call. = FALSE)
  for (a in jsonlite::fromJSON(path, simplifyDataFrame = FALSE)$args)
    p <- argparser::add_argument(p, a$flag, help = a$help %||% "",
           type = switch(a$type %||% "string", integer = "integer", number = "numeric", "character"))
  p
}

# On purpose, we don't try to enforce options or choices — left for later.

add_base_args  <- function(p) .inject(p, file.path(SCHEMA_DIR, "_base.json"))
add_stage_args <- function(p, schema) .inject(p, file.path(SCHEMA_DIR, paste0(schema, ".json")))
