#!/usr/bin/env Rscript
# sample_ena_runs.R — the ENA runs behind every sample, with sequencing depth and platform.
#
# Chain: publication sample code (data/open/samples.tsv) -> sequencing library name (data/quantification/
# sample_name_map.tsv) -> the lane file(s) the quantification used (data/quantification/lane_map.tsv, used == yes) ->
# the ENA run whose submitted file carries the same key. The key is the Novogene library prefix
# <library>_<EKRN id>_<flowcell>_<lane> for Arabidopsis, aspen and birch, and the facility library number plus its
# sample-sheet index (P#####_###_S#) for pine. A lane file built from all of its library's runs (a "_merge_" file, or
# the "TechRepComb_" spruce files that combine technical replicates) is matched on the library number alone, so it
# takes every run of that library: each spruce sample has two runs, and pine sample 120 was merged from two.
#
# Run metadata come from the ENA portal API (read_count = read pairs, instrument_model, submitted and FASTQ file
# names). For each run the first read of each FASTQ file (mate 1 and mate 2) is read from the first 64 kB of the file:
# the sequencer's instrument and flowcell fields from its header, and the read length from its sequence. ENA's
# instrument_model is the value submitted; the header is what the sequencer wrote. (ENA's base_count is not used: for
# some of these runs it counts one mate only.)
#
# Batches: runs are grouped by library naming scheme, and every run in a batch must share one instrument-ID prefix in
# its read headers (the script stops otherwise). Batch labels are those of the Methods (sequencing provider and
# platform).
#
# Inputs (deposit v4 or later, unpacked into ./data):
#   data/open/samples.tsv, data/quantification/sample_name_map.tsv, data/quantification/lane_map.tsv
# Outputs:
#   results/ena_run_metadata.tsv   one row per run: run, study, sample_id, species, library, read_pairs, read_length,
#                                  ena_instrument_model, header_instrument, header_flowcell, batch
#   results/sample_ena_runs.tsv    one row per sample: sample_id, species, ena_runs (";"-separated), read_pairs, batch,
#                                  runs_of_same_library_not_used (in ENA, but not in the quantification)
#   results/sequencing_depth.tsv   per batch: samples, runs, median read pairs per sample, read length(s),
#                                  ENA instrument_model value(s), header instrument prefix(es)
# Usage: Rscript scripts/deposit/sample_ena_runs.R

suppressMessages(library(curl))
STUDIES <- c("PRJEB104117", "PRJEB104118", "PRJEB104119", "PRJEB104120", "PRJEB26918")   # Data availability
API <- "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS <- "run_accession,study_accession,read_count,instrument_model,submitted_ftp,fastq_ftp"
BATCH_LABEL <- c(novogene = "Novogene NovaSeq 6000", P21452 = "NGI NovaSeq 6000", P1406 = "HiSeq 2500")

rd <- function(p) read.delim(p, colClasses = "character", check.names = FALSE, quote = "")
get_text <- function(url) rawToChar(curl_fetch_memory(url)$content)

ena_runs <- function() {
  do.call(rbind, lapply(STUDIES, function(st) {
    q <- sprintf("%s?result=read_run&query=study_accession%%3D%s&fields=%s&format=tsv&limit=0", API, st, FIELDS)
    read.delim(text = get_text(q), colClasses = "character", quote = "")
  }))
}

first_record <- function(url) {                       # header and sequence of the first read, from the first 64 kB
  h <- new_handle(); handle_setheaders(h, Range = "bytes=0-65535")
  raw <- curl_fetch_memory(paste0("https://", url), handle = h)$content
  con <- gzcon(rawConnection(raw)); on.exit(close(con))
  readLines(con, n = 2, warn = FALSE)
}

sequencer_fields <- function(fastq_ftp) {             # instrument, flowcell, read lengths ("150+150")
  recs <- lapply(strsplit(fastq_ftp, ";")[[1]], first_record)
  f <- lapply(recs, function(r) strsplit(strsplit(r[1], " ")[[1]][2], ":")[[1]])
  inst <- unique(vapply(f, function(x) paste(x[1], x[3]), ""))
  stopifnot(length(inst) == 1)
  x <- f[[1]]
  list(instrument = x[1], flowcell = x[3], read_length = paste(vapply(recs, function(r) nchar(r[2]), 0L), collapse = "+"))
}

key_of <- function(name, library_level = FALSE) {
  m <- regmatches(name, regexec("(P[0-9]{4,5}_[0-9]{3})(_S[0-9]+)?", name))[[1]]
  if (length(m)) return(if (library_level || m[3] == "") m[2] else paste0(m[2], m[3]))
  m <- regmatches(name, regexec("^(.+?_EKRN[0-9]+-[A-Za-z0-9_]+?_[A-Za-z0-9_]+?_L[0-9])", name))[[1]]
  if (length(m)) m[2] else NA_character_
}
combines_runs <- function(f) grepl("_merge_", f, fixed = TRUE) || startsWith(f, "TechRepComb_")
batch_of <- function(k) if (grepl("_EKRN", k, fixed = TRUE)) "novogene" else sub("_.*", "", k)
submitted_name <- function(r) basename(strsplit(r, ";")[[1]][1])

samples <- rd("data/open/samples.tsv")
nm <- rd("data/quantification/sample_name_map.tsv")
name <- setNames(nm$sequencing_name, paste(nm$genotype, nm$deposited_sample_id))
lm <- rd("data/quantification/lane_map.tsv"); lm <- lm[lm$used == "yes", ]

runs <- ena_runs()
runs$file <- vapply(runs$submitted_ftp, submitted_name, "")
runs$key <- vapply(runs$file, key_of, "")
runs$lib <- vapply(runs$file, key_of, "", library_level = TRUE)

run_rows <- list(); sample_rows <- list(); used <- character(0)
for (i in seq_len(nrow(samples))) {
  sp <- samples$species[i]; sid <- samples$sample_id[i]; lib <- name[[paste(sp, sid)]]
  lf <- lm$quantification_file[lm$genotype == sp & lm$deposited_sample_id == lib]
  whole <- vapply(lf, combines_runs, TRUE)
  keys <- unique(data.frame(k = mapply(key_of, lf, whole), whole = whole))
  if (!nrow(keys) || anyNA(keys$k)) stop(sp, " ", sid, ": no library key for lane files ", paste(lf, collapse = ", "))
  rr <- do.call(rbind, lapply(seq_len(nrow(keys)), function(j)
    runs[if (keys$whole[j]) runs$lib %in% keys$k[j] else runs$key %in% keys$k[j], ]))
  if (is.null(rr) || !nrow(rr)) stop(sp, " ", sid, ": no ENA run for ", paste(keys$k, collapse = ", "))
  b <- unique(vapply(keys$k, batch_of, "")); stopifnot(length(b) == 1)
  left_out <- sort(setdiff(runs$run_accession[runs$lib %in% unique(rr$lib)], rr$run_accession))
  used <- c(used, rr$run_accession)
  for (j in seq_len(nrow(rr))) {
    s <- sequencer_fields(rr$fastq_ftp[j])
    run_rows[[length(run_rows) + 1]] <- data.frame(
      run = rr$run_accession[j], study = rr$study_accession[j], sample_id = sid, species = sp, library = lib,
      read_pairs = as.numeric(rr$read_count[j]), read_length = s$read_length,
      ena_instrument_model = rr$instrument_model[j], header_instrument = s$instrument, header_flowcell = s$flowcell,
      batch = BATCH_LABEL[[b]])
  }
  sample_rows[[length(sample_rows) + 1]] <- data.frame(
    sample_id = sid, species = sp, ena_runs = paste(sort(rr$run_accession), collapse = ";"),
    read_pairs = sum(as.numeric(rr$read_count)), batch = BATCH_LABEL[[b]],
    runs_of_same_library_not_used = paste(left_out, collapse = ";"))
}
run_rows <- do.call(rbind, run_rows); sample_rows <- do.call(rbind, sample_rows)

prefix <- sub("^([A-Za-z-]+).*", "\\1", run_rows$header_instrument)
pref <- tapply(prefix, run_rows$batch, function(x) sort(unique(x)))
bad <- names(pref)[lengths(pref) != 1]
if (length(bad)) stop("batch ", bad[1], ": read headers carry more than one instrument prefix")

dir.create("results", showWarnings = FALSE)
wr <- function(x, p) write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
wr(run_rows, "results/ena_run_metadata.tsv")
wr(sample_rows, "results/sample_ena_runs.tsv")
depth <- do.call(rbind, lapply(sort(unique(sample_rows$batch)), function(b) {
  ss <- sample_rows[sample_rows$batch == b, ]; rr <- run_rows[run_rows$batch == b, ]
  data.frame(batch = b, samples = nrow(ss), runs = nrow(rr), median_read_pairs_per_sample = median(ss$read_pairs),
             read_lengths = paste(sort(unique(rr$read_length)), collapse = ";"),
             ena_instrument_models = paste(sort(unique(rr$ena_instrument_model)), collapse = ";"),
             header_instrument_prefix = paste(pref[[b]], collapse = ";"))
}))
wr(depth, "results/sequencing_depth.tsv")
unused <- sort(setdiff(runs$run_accession, used))
cat(sprintf("%d samples, %d runs used, %d runs in these studies not used\n", nrow(sample_rows), nrow(run_rows),
            length(unused)))
print(depth, row.names = FALSE)
