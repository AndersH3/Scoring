# NetCDF-4 storage for the Scoring project.
# All scientific data, checkpoints and numerical metadata are persisted in .nc.
# Source code, rendered reports and figures are NOT scientific data containers.
# Requires: ncdf4, jsonlite. No CSV/RDS output is permitted.
#
# Schema scoring-netcdf-v1:
#   global scoring_schema and scoring_tree_json attributes
#   one typed variable v000001, ... per atomic leaf
#   optional v000001_missing byte mask (1 = missing)
#   NC_CHAR strings with independent UTF-8 string-length dimension
#   numeric array dimensions are recorded in scoring_tree_json
#   tables are encoded as ordered, named vectors; nested lists are supported.
# Python readers need to reverse dimension order for multidimensional values,
# because ncdf4 uses R/Fortran dimension order on top of NetCDF's C order.

scoring_nc_check <- function() {
  for (p in c("ncdf4", "jsonlite")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("NetCDF dependency missing: install.packages('%s')", p),
           call. = FALSE)
  }
}

scoring_nc_path <- function(path) {
  if (!is.character(path) || length(path) != 1L ||
      !grepl("\\.nc$", path, ignore.case = TRUE))
    stop("Persisted research data must have a .nc extension", call. = FALSE)
  path
}

write_scoring_netcdf <- function(object, path, metadata = list()) {
  scoring_nc_check()
  scoring_nc_path(path)
  leaves <- list()
  build <- function(x) {
    if (is.data.frame(x) || (is.list(x) && !is.object(x))) {
      kind <- if (is.data.frame(x)) "data.frame" else "list"
      return(list(kind = kind, names = as.list(names(x) %||% rep("", length(x))),
                  entries = unname(lapply(x, build)),
                  nrow = if (is.data.frame(x)) nrow(x) else NULL))
    }
    if (is.null(x)) return(list(kind = "null"))
    if (is.factor(x)) x <- as.character(x)
    if (inherits(x, "Date") || inherits(x, "POSIXt")) x <- as.character(x)
    if (!(is.integer(x) || is.numeric(x) || is.logical(x) || is.character(x)))
      stop("Unsupported NetCDF value type: ", paste(class(x), collapse = "/"),
           call. = FALSE)
    kind <- if (is.character(x)) "string" else if (is.logical(x)) "logical" else
            if (is.integer(x)) "integer" else "double"
    id <- sprintf("v%06d", length(leaves) + 1L)
    shape <- dim(x)
    if (is.null(shape)) shape <- as.integer(length(x))
    has_na <- anyNA(x)
    leaves[[length(leaves) + 1L]] <<-
      list(id = id, value = x, shape = as.integer(shape), kind = kind, has_na = has_na)
    list(kind = "leaf", variable = id, type = kind,
         dimensions = as.list(as.integer(shape)), has_na = has_na,
         element_names = if (is.null(names(x))) NULL else as.list(names(x)))
  }
  "%||%" <- function(a, b) if (is.null(a)) b else a
  tree <- build(object)
  defs <- list()
  for (leaf in leaves) {
    if (length(leaf$value) == 0L)
      stop("Empty atomic leaves must be represented by NULL", call. = FALSE)
    dims <- lapply(seq_along(leaf$shape), function(k)
      ncdf4::ncdim_def(sprintf("dim_%s_%d", leaf$id, k), "",
                      seq_len(leaf$shape[[k]]), create_dimvar = FALSE))
    leaf$dimensions <- dims
    prec <- switch(leaf$kind, string = "char", integer = "integer",
                   logical = "byte", double = "double")
    if (identical(leaf$kind, "string")) {
      # nchar(type='bytes') counts UTF-8 octets, NOT Unicode code points.
      s <- enc2utf8(replace(leaf$value, is.na(leaf$value), ""))
      width <- max(1L, nchar(s, type = "bytes"))
      strlen_dim <- ncdf4::ncdim_def(paste0("strlen_", leaf$id), "",
                                    seq_len(width), create_dimvar = FALSE)
      dims <- c(list(strlen_dim), dims)
      leaf$value <- s
    }
    defs[[length(defs) + 1L]] <- ncdf4::ncvar_def(
      leaf$id, "", dims, prec = prec, missval = NA,
      compression = 4, shuffle = TRUE)
    if (leaf$has_na) {
      defs[[length(defs) + 1L]] <- ncdf4::ncvar_def(
        paste0(leaf$id, "_missing"), "", leaf$dimensions,
        prec = "byte", missval = NA, compression = 4)
    }
    leaves[[as.integer(substring(leaf$id, 2L))]] <- leaf
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".scoring-", tmpdir = dirname(path), fileext = ".nc")
  on.exit(unlink(tmp), add = TRUE)
  nc <- ncdf4::nc_create(tmp, defs, force_v4 = TRUE)
  closed <- FALSE
  on.exit(if (!closed) ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncatt_put(nc, 0, "scoring_schema", "scoring-netcdf-v1")
  ncdf4::ncatt_put(nc, 0, "scoring_tree_json",
                   as.character(jsonlite::toJSON(tree, auto_unbox = TRUE,
                                                  null = "null", digits = NA)))
  ncdf4::ncatt_put(nc, 0, "scoring_metadata_json",
                   as.character(jsonlite::toJSON(metadata, auto_unbox = TRUE,
                                                  null = "null", digits = NA)))
  for (leaf in leaves) {
    v <- leaf$value
    if (leaf$kind != "string") {
      v[is.na(v)] <- 0
      if (leaf$kind == "logical") v <- as.integer(v)
    }
    ncdf4::ncvar_put(nc, leaf$id, v)
    if (leaf$has_na)
      ncdf4::ncvar_put(nc, paste0(leaf$id, "_missing"),
                       as.integer(is.na(leaves[[as.integer(substring(leaf$id,2L))]]$value)))
  }
  ncdf4::nc_sync(nc)
  ncdf4::nc_close(nc)
  closed <- TRUE
  # A failed rename must never be described as a completed save.
  if (!file.rename(tmp, path))
    stop("Failed to atomically replace NetCDF output: ", path, call. = FALSE)
  invisible(path)
}

read_scoring_netcdf <- function(path) {
  scoring_nc_check()
  scoring_nc_path(path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  schema <- ncdf4::ncatt_get(nc, 0, "scoring_schema")
  if (!schema$hasatt || !identical(schema$value, "scoring-netcdf-v1"))
    stop("File is not a scoring-netcdf-v1 NetCDF file", call. = FALSE)
  tree <- jsonlite::fromJSON(ncdf4::ncatt_get(nc, 0, "scoring_tree_json")$value,
                             simplifyVector = FALSE)
  rebuild <- function(node) {
    kind <- node$kind
    if (kind == "null") return(NULL)
    if (kind == "leaf") {
      val <- ncdf4::ncvar_get(nc, node$variable, collapse_degen = TRUE)
      shape <- as.integer(unlist(node$dimensions))
      if (node$type == "string") val <- as.character(val)
      if (node$type == "integer") val <- as.integer(val)
      if (node$type == "logical") val <- as.logical(val)
      if (node$type == "double") val <- as.double(val)
      if (isTRUE(node$has_na)) {
        missing <- as.logical(ncdf4::ncvar_get(
          nc, paste0(node$variable, "_missing"), collapse_degen = TRUE))
        val[missing] <- NA
      }
      if (length(shape) > 1L) dim(val) <- shape
      if (!is.null(node$element_names)) names(val) <- unlist(node$element_names)
      return(val)
    }
    entries <- lapply(node$entries, rebuild)
    nm <- unlist(node$names)
    if (length(nm) == length(entries)) names(entries) <- nm
    if (kind == "data.frame") {
      result <- as.data.frame(entries, optional = TRUE,
                              stringsAsFactors = FALSE, check.names = FALSE)
      names(result) <- nm
      return(result)
    }
    if (kind == "list") return(entries)
    stop("Unknown NetCDF schema node: ", kind, call. = FALSE)
  }
  rebuild(tree)
}

write_scoring_table <- function(data, path, metadata = list()) {
  if (!is.data.frame(data)) stop("Expected a data.frame", call. = FALSE)
  write_scoring_netcdf(data, path, metadata)
}

read_scoring_table <- function(path) {
  result <- read_scoring_netcdf(path)
  if (!is.data.frame(result)) stop("Expected NetCDF table", call. = FALSE)
  result
}

# Reject non-NetCDF persistence, even if a caller supplies a .csv/.rds name.
# Legacy CSV may be imported in memory by a deliberate migration utility only.
