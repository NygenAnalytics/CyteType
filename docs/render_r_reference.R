#!/usr/bin/env Rscript
# Render the R function reference, docs/source/reference/r-api.md, from the
# help pages (man/*.Rd) of the CyteTypeR package.
#
# Usage, from the root of the CyteType repository:
#
#   Rscript docs/render_r_reference.R [path/to/CyteTypeR]
#
# The path defaults to ../CyteTypeR, a checkout of CyteTypeR next to this
# repository. Run roxygen2::roxygenise() in CyteTypeR first if its R sources
# changed. The script uses base R only (tools::parse_Rd) and installs nothing.
# It writes an H2 section for each group in GROUPS below and, under it, an H3
# section for each of the group's exported functions; add new exported
# functions to a group there.

GROUPS <- list(
  list(
    title = "Prepare and run",
    label = "r-api-prepare-and-run",
    intro = paste(
      "Prepare a Seurat object, submit the annotation job and wait for the",
      "results, then delete the uploaded files."
    ),
    functions = c("PrepareCyteTypeR", "CyteTypeR", "CleanUpArtifacts")
  ),
  list(
    title = "Results",
    label = "r-api-results",
    intro = paste(
      "Read the results stored in a Seurat object, or fetch the results of a",
      "job by its ID."
    ),
    functions = c("GetResults")
  ),
  list(
    title = "Sign-in and the cytetyper command",
    label = "r-api-sign-in",
    intro = paste(
      "Save an API key, open the dashboard or a job's report, and install the",
      "`cytetyper` terminal command. Each `cytetyper` command calls one of",
      "these functions (see the [Command-line reference](cli.md))."
    ),
    functions = c(
      "SetupCyteTypeR",
      "LoginCyteTypeR",
      "LogoutCyteTypeR",
      "OpenCyteTypeDashboard",
      "ViewCyteTypeJob",
      "InstallCyteTypeRCli"
    )
  ),
  list(
    title = "Errors",
    label = "r-api-errors",
    intro = paste(
      "The condition that CyteTypeR signals when the status or the results of",
      "a job cannot be fetched (see [Errors and error codes](errors.md))."
    ),
    functions = c("cytetype_api_error")
  )
)

DESCRIPTION_LINE <- paste(
  "description: The exported functions of the CyteTypeR R package, with their",
  "usage, arguments, return values and examples, generated from the package's",
  "help pages."
)

# Paths ---------------------------------------------------------------------

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath("docs/render_r_reference.R", mustWork = FALSE))
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

repository_root <- dirname(dirname(script_path()))
arguments <- commandArgs(trailingOnly = TRUE)
package_path <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  file.path(repository_root, "..", "CyteTypeR")
}
package_path <- normalizePath(package_path, mustWork = FALSE)
if (!file.exists(file.path(package_path, "DESCRIPTION")) ||
    !dir.exists(file.path(package_path, "man"))) {
  stop("No CyteTypeR package (DESCRIPTION and man/) at ", package_path, call. = FALSE)
}
output_path <- file.path(repository_root, "docs", "source", "reference", "r-api.md")

description <- read.dcf(file.path(package_path, "DESCRIPTION"), fields = c("Package", "Version"))
if (!identical(unname(description[1L, "Package"]), "CyteTypeR")) {
  stop(package_path, " is not the CyteTypeR package", call. = FALSE)
}
package_version <- unname(description[1L, "Version"])

# Exported functions from NAMESPACE, in either roxygen2 layout.
namespace_text <- paste(readLines(file.path(package_path, "NAMESPACE"), warn = FALSE), collapse = "\n")
export_blocks <- regmatches(namespace_text, gregexpr("(?m)^export\\(([^)]*)\\)", namespace_text, perl = TRUE))[[1L]]
exports <- unlist(lapply(export_blocks, function(block) {
  inner <- sub("^export\\(", "", sub("\\)$", "", block))
  trimws(gsub("[\"`]", "", strsplit(inner, ",", fixed = TRUE)[[1L]]))
}))
exports <- exports[nzchar(exports)]

# Rd helpers -----------------------------------------------------------------

tag <- function(node) {
  value <- attr(node, "Rd_tag")
  if (is.null(value)) "" else value
}

section <- function(rd, name) {
  matches <- Filter(function(node) identical(tag(node), name), rd)
  if (length(matches) == 0L) NULL else matches[[1L]]
}

sections <- function(rd, name) {
  Filter(function(node) identical(tag(node), name), rd)
}

anchor <- function(name) {
  paste0("r-api-", gsub("^-+|-+$", "", gsub("[^a-z0-9]+", "-", tolower(name))))
}

# Plain text of an Rd fragment, as written (for code, usage and examples).
plain_text <- function(node) {
  node_tag <- tag(node)
  if (node_tag %in% c("TEXT", "RCODE", "VERB")) {
    return(paste(as.character(node), collapse = ""))
  }
  if (node_tag == "COMMENT") {
    return("")
  }
  if (node_tag %in% c("\\dots", "\\ldots")) {
    return("...")
  }
  if (node_tag %in% c("\\dontshow", "\\testonly")) {
    return("")
  }
  if (node_tag == "\\method" || node_tag == "\\S3method") {
    return(paste0(plain_text(node[[1L]])))
  }
  if (node_tag == "\\R") {
    return("R")
  }
  if (is.list(node)) {
    return(paste(vapply(node, plain_text, character(1)), collapse = ""))
  }
  paste(as.character(node), collapse = "")
}

escape_text <- function(text) {
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  text <- gsub("([`*_\\[\\]<>$|~])", "\\\\\\1", text, perl = TRUE)
  gsub("(?m)^([ \t]*)#", "\\1\\\\#", text, perl = TRUE)
}

code_span <- function(text) {
  text <- gsub("\\s+", " ", text)
  fence <- if (grepl("`", text, fixed = TRUE)) "``" else "`"
  paste0(fence, text, fence)
}

link_target <- function(node) {
  option <- attr(node, "Rd_option")
  topic <- plain_text(node)
  if (!is.null(option)) {
    option <- plain_text(option)
    if (startsWith(option, "=")) {
      topic <- substring(option, 2L)
    } else {
      return(NULL)
    }
  }
  topic <- sub("\\(\\)$", "", topic)
  if (topic %in% exports) topic else NULL
}

# Inline Markdown of an Rd fragment. Newlines are kept so that paragraphs can
# be split later; lists are flattened when `flatten_lists` is TRUE.
inline <- function(node, flatten_lists = FALSE) {
  node_tag <- tag(node)
  recurse <- function(x) {
    paste(vapply(x, inline, character(1), flatten_lists = flatten_lists), collapse = "")
  }
  switch(
    node_tag,
    "TEXT" = escape_text(paste(as.character(node), collapse = "")),
    "RCODE" = ,
    "VERB" = escape_text(paste(as.character(node), collapse = "")),
    "COMMENT" = "",
    "\\code" = ,
    "\\verb" = ,
    "\\samp" = ,
    "\\file" = ,
    "\\env" = ,
    "\\option" = ,
    "\\command" = ,
    "\\kbd" = {
      links <- Filter(function(child) identical(tag(child), "\\link"), node)
      text <- plain_text(node)
      target <- if (length(links) == 1L) link_target(links[[1L]]) else NULL
      if (is.null(target)) {
        code_span(text)
      } else {
        paste0("[", code_span(text), "](#", anchor(target), ")")
      }
    },
    "\\link" = {
      target <- link_target(node)
      text <- plain_text(node)
      if (is.null(target)) escape_text(text) else paste0("[", escape_text(text), "](#", anchor(target), ")")
    },
    "\\url" = {
      url <- trimws(plain_text(node))
      paste0("[", escape_text(url), "](", url, ")")
    },
    "\\href" = paste0("[", recurse(node[[2L]]), "](", trimws(plain_text(node[[1L]])), ")"),
    "\\email" = {
      email <- trimws(plain_text(node))
      paste0("[", escape_text(email), "](mailto:", email, ")")
    },
    "\\strong" = ,
    "\\bold" = paste0("**", recurse(node), "**"),
    "\\emph" = ,
    "\\var" = ,
    "\\dfn" = paste0("*", recurse(node), "*"),
    "\\sQuote" = paste0("'", recurse(node), "'"),
    "\\dQuote" = paste0("\"", recurse(node), "\""),
    "\\pkg" = ,
    "\\acronym" = ,
    "\\abbr" = ,
    "\\cite" = recurse(node),
    "\\dots" = ,
    "\\ldots" = "...",
    "\\R" = "R",
    "\\cr" = "\n",
    "\\enc" = recurse(node[[1L]]),
    "\\eqn" = code_span(plain_text(node[[length(node)]])),
    "\\if" = {
      format <- plain_text(node[[1L]])
      if (grepl("html|text", format)) recurse(node[[2L]]) else ""
    },
    "\\ifelse" = {
      format <- plain_text(node[[1L]])
      if (grepl("html|text", format)) recurse(node[[2L]]) else recurse(node[[3L]])
    },
    "\\out" = plain_text(node),
    "\\Sexpr" = ,
    "\\figure" = ,
    "\\special" = "",
    "\\itemize" = ,
    "\\enumerate" = ,
    "\\describe" = {
      items <- list_items(node)
      if (flatten_lists) {
        paste0(" ", paste(items, collapse = "; "), " ")
      } else {
        stop("Block element inside inline content: ", node_tag, call. = FALSE)
      }
    },
    if (is.list(node)) recurse(node) else escape_text(paste(as.character(node), collapse = ""))
  )
}

squash <- function(text) {
  trimws(gsub("[ \t]*\n[ \t]*", " ", gsub("[ \t]+", " ", text)))
}

# The items of \itemize, \enumerate or \describe as inline Markdown.
list_items <- function(node) {
  node_tag <- tag(node)
  if (node_tag == "\\describe") {
    items <- Filter(function(child) identical(tag(child), "\\item"), node)
    return(vapply(items, function(item) {
      paste0("**", squash(inline(item[[1L]], TRUE)), "**: ", squash(inline(item[[2L]], TRUE)))
    }, character(1)))
  }
  items <- list()
  current <- NULL
  for (child in node) {
    if (identical(tag(child), "\\item")) {
      if (!is.null(current)) items[[length(items) + 1L]] <- current
      current <- list()
    } else if (!is.null(current)) {
      current[[length(current) + 1L]] <- child
    }
  }
  if (!is.null(current)) items[[length(items) + 1L]] <- current
  vapply(items, function(item) {
    squash(paste(vapply(item, inline, character(1), flatten_lists = TRUE), collapse = ""))
  }, character(1))
}

# Block Markdown (paragraphs, lists and code blocks) of an Rd section.
blocks <- function(node) {
  output <- character()
  buffer <- ""
  flush <- function() {
    paragraphs <- strsplit(buffer, "\n[ \t]*\n", perl = TRUE)[[1L]]
    paragraphs <- vapply(paragraphs, squash, character(1), USE.NAMES = FALSE)
    output <<- c(output, paragraphs[nzchar(paragraphs)])
    buffer <<- ""
  }
  for (child in node) {
    child_tag <- tag(child)
    if (child_tag %in% c("\\itemize", "\\enumerate", "\\describe")) {
      flush()
      items <- list_items(child)
      marker <- if (child_tag == "\\enumerate") paste0(seq_along(items), ". ") else rep("- ", length(items))
      output <- c(output, paste(paste0(marker, items), collapse = "\n"))
    } else if (child_tag %in% c("\\preformatted", "\\deqn")) {
      flush()
      output <- c(output, paste0("```text\n", trim_code(plain_text(child)), "\n```"))
    } else if (child_tag == "\\item") {
      # \item{name}{description} entries of \value.
      flush()
      output <- c(output, paste0("- ", code_span(plain_text(child[[1L]])), ": ", squash(inline(child[[2L]], TRUE))))
    } else {
      buffer <- paste0(buffer, inline(child))
    }
  }
  flush()
  paste(output, collapse = "\n\n")
}

trim_code <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  lines <- sub("[ \t]+$", "", lines)
  while (length(lines) > 0L && !nzchar(lines[[1L]])) lines <- lines[-1L]
  while (length(lines) > 0L && !nzchar(lines[[length(lines)]])) lines <- lines[-length(lines)]
  indents <- nchar(lines[nzchar(lines)]) - nchar(sub("^ +", "", lines[nzchar(lines)]))
  if (length(indents) > 0L && min(indents) > 0L) {
    lines <- substring(lines, min(indents) + 1L)
  }
  paste(lines, collapse = "\n")
}

table_cell <- function(text) {
  gsub("\n", " ", text, fixed = TRUE)
}

# One function -----------------------------------------------------------------

render_function <- function(name, rd) {
  title <- squash(inline(section(rd, "\\title"), TRUE))
  description_text <- blocks(section(rd, "\\description"))
  lines <- c(
    paste0("(", anchor(name), ")="),
    paste0("### `", name, "()`"),
    ""
  )
  if (!identical(squash(gsub("[.]$", "", title)), squash(gsub("[.]$", "", description_text)))) {
    lines <- c(lines, paste0("*", title, "*"), "")
  }
  lines <- c(lines, description_text, "")

  usage <- section(rd, "\\usage")
  if (!is.null(usage)) {
    lines <- c(lines, "**Usage**", "", "```r", trim_code(plain_text(usage)), "```", "")
  }

  arguments <- section(rd, "\\arguments")
  if (!is.null(arguments)) {
    items <- Filter(function(child) identical(tag(child), "\\item"), arguments)
    rows <- vapply(items, function(item) {
      argument_names <- trimws(strsplit(plain_text(item[[1L]]), ",", fixed = TRUE)[[1L]])
      paste0(
        "| ", paste(vapply(argument_names, code_span, character(1)), collapse = ", "),
        " | ", table_cell(squash(inline(item[[2L]], TRUE))), " |"
      )
    }, character(1))
    lines <- c(lines, "**Arguments**", "", "| Argument | Description |", "| --- | --- |", rows, "")
  }

  value <- section(rd, "\\value")
  if (!is.null(value)) {
    lines <- c(lines, "**Value**", "", blocks(value), "")
  }

  details <- section(rd, "\\details")
  if (!is.null(details)) {
    lines <- c(lines, "**Details**", "", blocks(details), "")
  }

  for (custom in sections(rd, "\\section")) {
    lines <- c(lines, paste0("**", squash(inline(custom[[1L]], TRUE)), "**"), "", blocks(custom[[2L]]), "")
  }

  examples <- section(rd, "\\examples")
  if (!is.null(examples)) {
    lines <- c(lines, "**Examples**", "", "```r", trim_code(plain_text(examples)), "```", "")
  }

  see_also <- section(rd, "\\seealso")
  if (!is.null(see_also)) {
    lines <- c(lines, "**See also**", "", blocks(see_also), "")
  }
  lines
}

# Page -----------------------------------------------------------------------

rd_files <- sort(list.files(file.path(package_path, "man"), pattern = "[.]Rd$", full.names = TRUE))
pages <- lapply(rd_files, tools::parse_Rd, permissive = FALSE)
documented <- list()
for (rd in pages) {
  doc_type <- section(rd, "\\docType")
  if (!is.null(doc_type) && identical(trimws(plain_text(doc_type)), "package")) {
    next
  }
  aliases <- vapply(sections(rd, "\\alias"), function(node) trimws(plain_text(node)), character(1))
  for (alias in intersect(aliases, exports)) {
    documented[[alias]] <- rd
  }
}

grouped <- unlist(lapply(GROUPS, `[[`, "functions"))
missing_docs <- setdiff(grouped, names(documented))
if (length(missing_docs) > 0L) {
  stop("No help page for: ", paste(missing_docs, collapse = ", "), call. = FALSE)
}
ungrouped <- setdiff(names(documented), grouped)
if (length(ungrouped) > 0L) {
  warning(
    "Exported functions without a group are listed under 'Other functions': ",
    paste(ungrouped, collapse = ", "),
    call. = FALSE
  )
  GROUPS[[length(GROUPS) + 1L]] <- list(
    title = "Other functions",
    label = "r-api-other",
    intro = "Other exported functions of CyteTypeR.",
    functions = sort(ungrouped)
  )
}

summary_rows <- unlist(lapply(GROUPS, function(group) {
  vapply(group$functions, function(name) {
    title <- squash(inline(section(documented[[name]], "\\title"), TRUE))
    paste0("| [`", name, "()`](#", anchor(name), ") | ", group$title, " | ", table_cell(title), " |")
  }, character(1))
}))

page <- c(
  "---",
  DESCRIPTION_LINE,
  "---",
  "",
  "(r-api)=",
  "# R function reference",
  "",
  paste0(
    "This page is generated from the help pages of CyteTypeR ", package_version,
    ". To regenerate it, run `Rscript docs/render_r_reference.R ../CyteTypeR` from",
    " the root of the CyteType repository. In R, the same help is available with",
    " `?` and the function name, for example `?PrepareCyteTypeR`."
  ),
  "",
  "| Function | Section | Purpose |",
  "| --- | --- | --- |",
  summary_rows,
  ""
)
for (group in GROUPS) {
  page <- c(
    page,
    paste0("(", group$label, ")="),
    paste0("## ", group$title),
    "",
    group$intro,
    ""
  )
  for (name in group$functions) {
    page <- c(page, render_function(name, documented[[name]]))
  }
}

page <- sub("[ \t]+$", "", page)
text <- paste(page, collapse = "\n")
text <- gsub("\n{3,}", "\n\n", text)
if (grepl("<!--", text, fixed = TRUE)) {
  stop("The rendered page contains an HTML comment", call. = FALSE)
}
writeLines(enc2utf8(text), output_path, useBytes = TRUE)
message("Wrote ", output_path, " from CyteTypeR ", package_version, " (", length(documented), " functions).")
