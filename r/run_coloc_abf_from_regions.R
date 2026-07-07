# Exploratory COLOC runner using regional TSVs produced by the stronger-genetics stack.
# This script is designed to be *runnable* on Windows with base R + coloc installed.
#
# IMPORTANT:
# - Requires allele harmonization. We do a best-effort harmonization using:
#     GWAS: A_EFF / A_NONEFF
#     eQTL: variant string "chrN_pos_REF_ALT" (assume beta is ALT effect)
# - Case/control fraction is not available in Bakker regional TSVs; we treat GWAS as
#   quantitative with N = Neff for an exploratory PP4 screen. For publication-grade
#   inference, provide s (case fraction) and use type="cc" or use SMR/HEIDI / SuSiE-coloc.
#
# Usage:
#   Rscript r/run_coloc_abf_from_regions.R --locus-id MARCO_rs2422391
#   Rscript r/run_coloc_abf_from_regions.R --tier1-tsv results/phases_g/tier1_science_supplement.tsv --top-n 10
#   Rscript r/run_coloc_abf_from_regions.R --locus-ids FN3K_rs1046896,TRPC4AP_rs73903009 --chrpos-fallback --merge-summary
# Flags:
#   --chrpos-fallback   If rsID merge + strict alleles yields 0 SNPs, merge on CHR:BP (same norm as GWAS/eQTL TSVs).
#   --relaxed-alleles   If overlap is still tiny/empty, keep rsID (or chr:pos) overlaps with finite betas/SEs only (screening).
#   --merge-summary     Merge new rows into existing coloc_summary.tsv (replace same locus_id).
#   --candidates-tsv    Table with locus_id + blood_eqtl_dataset_id (default: tier1_coloc_candidates.tsv).
#
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default=NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

has_flag <- function(flag) any(args == flag)

# Rscript does not always populate sys.frame(1)$ofile on Windows.
# Use the current working directory as the project root (this script is invoked
# from the repo root in our pipeline commands).
proj_dir <- normalizePath(getwd(), winslash="\\", mustWork=TRUE)

out_dir <- file.path(proj_dir, "results", "phases_g", "stronger_genetics", "coloc")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("coloc", quietly = TRUE)) {
  stop("R package `coloc` not installed. Install with install.packages('coloc').")
}

parse_eqtl_variant <- function(v) {
  # v like "chr2_119214898_G_T" or indels "chr2_119234176_C_CT" / multi-part ALT
  parts <- strsplit(as.character(v), "_", fixed=TRUE)[[1]]
  if (length(parts) < 4) return(list(chr=NA, pos=NA, ref=NA, alt=NA))
  chr <- gsub("^chr", "", parts[1], ignore.case=TRUE)
  pos <- suppressWarnings(as.integer(parts[2]))
  ref <- parts[3]
  alt <- if (length(parts) == 4) parts[4] else paste(parts[4:length(parts)], collapse="_")
  list(chr=chr, pos=pos, ref=ref, alt=alt)
}

norm_chr_lab <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- gsub("^CHR", "", s, ignore.case = TRUE)
  out <- s
  out[s %in% c("23", "X")] <- "X"
  mask_num <- !(s %in% c("23", "X"))
  if (any(mask_num)) {
    suppressWarnings(n <- as.integer(s[mask_num]))
    good <- !is.na(n)
    wh <- which(mask_num)
    out[wh[good]] <- as.character(n[good])
  }
  out
}

load_region_pair <- function(locus_id) {
  gwas_path <- file.path(proj_dir, "results", "phases_g", "stronger_genetics", paste0("gwas_region_", locus_id, ".tsv"))
  eqtl_path <- file.path(proj_dir, "results", "phases_g", "stronger_genetics", paste0("eqtl_region_", locus_id, ".tsv"))
  if (!file.exists(gwas_path)) stop(paste("Missing", gwas_path))
  if (!file.exists(eqtl_path)) stop(paste("Missing", eqtl_path))
  if (file.info(eqtl_path)$size < 10) stop(paste("Empty eQTL file:", eqtl_path))
  gwas <- read.table(gwas_path, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="")
  eqtl <- read.table(eqtl_path, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="")
  list(gwas=gwas, eqtl=eqtl, gwas_path=gwas_path, eqtl_path=eqtl_path)
}

dedup_rs_merge <- function(m) {
  if (nrow(m) == 0) return(m)
  if ("pvalue" %in% colnames(m)) {
    o <- order(m$SNP, m$pvalue, na.last=TRUE)
    m <- m[o, , drop=FALSE]
    m <- m[!duplicated(m$SNP), , drop=FALSE]
  } else {
    m <- m[!duplicated(m$SNP), , drop=FALSE]
  }
  m
}

filter_eqtl_dataset <- function(eqtl, dataset_id) {
  if (!is.null(dataset_id) && nzchar(as.character(dataset_id)[1]) && "dataset_id" %in% colnames(eqtl)) {
    sub <- eqtl[eqtl$dataset_id == dataset_id, , drop=FALSE]
    if (nrow(sub) > 0) return(sub)
  }
  eqtl
}

attach_eqtl_alleles <- function(eqtl) {
  parsed <- lapply(eqtl$variant, parse_eqtl_variant)
  eqtl$eqtl_ref <- vapply(parsed, function(x) x$ref, character(1))
  eqtl$eqtl_alt <- vapply(parsed, function(x) x$alt, character(1))
  eqtl
}

allele_align <- function(m, strict) {
  m$gwas_beta_aligned <- m$BETA
  flip <- (!is.na(m$eqtl_alt)) & (toupper(m$A_EFF) == toupper(m$eqtl_ref)) & (toupper(m$A_NONEFF) == toupper(m$eqtl_alt))
  ok   <- (!is.na(m$eqtl_alt)) & (toupper(m$A_EFF) == toupper(m$eqtl_alt)) & (toupper(m$A_NONEFF) == toupper(m$eqtl_ref))
  m$gwas_beta_aligned[flip] <- -m$BETA[flip]
  if (strict) {
    m <- m[ok | flip, , drop=FALSE]
  } else {
    keep <- is.finite(m$BETA) & is.finite(m$SE) & is.finite(m$beta) & is.finite(m$se)
    m <- m[keep, , drop=FALSE]
    if (nrow(m) > 0) {
      warning(paste("relaxed-alleles: kept", nrow(m), "SNPs without strict REF/ALT match; exploratory only"))
    }
  }
  m
}

harmonize_and_coloc <- function(locus_id, dataset_id=NULL) {
  pair <- load_region_pair(locus_id)
  gwas <- pair$gwas
  eqtl0 <- pair$eqtl
  eqtl <- filter_eqtl_dataset(eqtl0, dataset_id)

  if (nrow(eqtl) < 50) {
    warning(paste("Low eQTL rows for locus", locus_id, "after dataset filter"))
  }

  if (!all(c("SNP","BP","BETA","SE","A_EFF","A_NONEFF","CHR") %in% colnames(gwas))) {
    stop("GWAS region TSV missing required columns (SNP, CHR, BP, BETA, SE, A_EFF, A_NONEFF)")
  }
  if (!all(c("rsid","position","beta","se","variant","chromosome") %in% colnames(eqtl))) {
    stop("eQTL region TSV missing required columns (rsid, chromosome, position, beta, se, variant)")
  }

  eqtl <- attach_eqtl_alleles(eqtl)

  m_rs <- merge(gwas, eqtl, by.x="SNP", by.y="rsid")
  m_rs <- dedup_rs_merge(m_rs)

  m <- allele_align(m_rs, strict=TRUE)

  if (nrow(m) == 0 && has_flag("--chrpos-fallback")) {
    g1 <- gwas
    e1 <- attach_eqtl_alleles(filter_eqtl_dataset(eqtl0, dataset_id))
    g1$cp <- paste(norm_chr_lab(g1$CHR), g1$BP, sep=":")
    e1$cp <- paste(norm_chr_lab(e1$chromosome), e1$position, sep=":")
    m_cp <- merge(g1, e1, by="cp", suffixes=c(".gw", ".eq"))
    if ("SNP.gw" %in% colnames(m_cp)) {
      m_cp$SNP <- m_cp$SNP.gw
    } else if (!"SNP" %in% colnames(m_cp)) {
      m_cp$SNP <- m_cp$rsid
    }
    if ("BP.gw" %in% colnames(m_cp)) m_cp$BP <- m_cp$BP.gw
    if ("BETA.gw" %in% colnames(m_cp)) m_cp$BETA <- m_cp$BETA.gw
    if ("SE.gw" %in% colnames(m_cp)) m_cp$SE <- m_cp$SE.gw
    if ("A_EFF.gw" %in% colnames(m_cp)) m_cp$A_EFF <- m_cp$A_EFF.gw
    if ("A_NONEFF.gw" %in% colnames(m_cp)) m_cp$A_NONEFF <- m_cp$A_NONEFF.gw
    if ("beta.eq" %in% colnames(m_cp)) m_cp$beta <- m_cp$beta.eq
    if ("se.eq" %in% colnames(m_cp)) m_cp$se <- m_cp$se.eq
    if ("position.eq" %in% colnames(m_cp)) m_cp$position <- m_cp$position.eq
    m_cp <- dedup_rs_merge(m_cp)
    m <- allele_align(m_cp, strict=TRUE)
    if (has_flag("--debug")) cat("[DEBUG] chrpos-fallback merged SNPs:", nrow(m), "\n")
  }

  if ((nrow(m) == 0 || nrow(m) < 15) && has_flag("--relaxed-alleles")) {
    m <- allele_align(m_rs, strict=FALSE)
    if (nrow(m) == 0 && has_flag("--chrpos-fallback")) {
      g1 <- gwas
      e1 <- attach_eqtl_alleles(filter_eqtl_dataset(eqtl0, dataset_id))
      g1$cp <- paste(norm_chr_lab(g1$CHR), g1$BP, sep=":")
      e1$cp <- paste(norm_chr_lab(e1$chromosome), e1$position, sep=":")
      m_cp <- merge(g1, e1, by="cp", suffixes=c(".gw", ".eq"))
      if ("SNP.gw" %in% colnames(m_cp)) m_cp$SNP <- m_cp$SNP.gw
      if ("BP.gw" %in% colnames(m_cp)) m_cp$BP <- m_cp$BP.gw
      if ("BETA.gw" %in% colnames(m_cp)) m_cp$BETA <- m_cp$BETA.gw
      if ("SE.gw" %in% colnames(m_cp)) m_cp$SE <- m_cp$SE.gw
      if ("A_EFF.gw" %in% colnames(m_cp)) m_cp$A_EFF <- m_cp$A_EFF.gw
      if ("A_NONEFF.gw" %in% colnames(m_cp)) m_cp$A_NONEFF <- m_cp$A_NONEFF.gw
      if ("beta.eq" %in% colnames(m_cp)) m_cp$beta <- m_cp$beta.eq
      if ("se.eq" %in% colnames(m_cp)) m_cp$se <- m_cp$se.eq
      if ("position.eq" %in% colnames(m_cp)) m_cp$position <- m_cp$position.eq
      m_cp <- dedup_rs_merge(m_cp)
      m <- allele_align(m_cp, strict=FALSE)
    }
  }

  if (nrow(m) < 50) {
    warning(paste("After harmonization, overlap SNPs for locus", locus_id, ":", nrow(m)))
  }
  if (nrow(m) == 0) stop("No overlapping SNPs after SNP / harmonization alignment")

  # Build coloc inputs
  # For a *programmatic* first-pass run, we set sdY=1 for both traits.
  # This avoids requiring per-study sdY or accurate case/control fraction.
  # For publication-grade analysis, replace this with trait-specific sdY (or MAF+N for sdY estimation)
  # and use type="cc" for the IA GWAS if case fraction is known.
  d1 <- list(
    beta = m$gwas_beta_aligned,
    varbeta = (m$SE)^2,
    snp = m$SNP,
    position = m$BP,
    type = "quant",
    sdY = 1
  )
  d2 <- list(
    beta = m$beta,
    varbeta = (m$se)^2,
    snp = m$SNP,
    position = m$position,
    type = "quant",
    sdY = 1
  )

  # Ensure identical SNP set and order (coloc requirement)
  ord <- match(d1$snp, d2$snp)
  keep <- which(!is.na(ord))
  if (length(keep) < 50) {
    warning(paste("After SNP alignment, overlap SNPs for locus", locus_id, ":", length(keep)))
  }
  if (length(keep) == 0) stop("No overlapping SNPs after SNP alignment")
  ord <- ord[keep]
  keep_vec <- function(x, idx) {
    if (is.null(x)) return(NULL)
    x[idx]
  }
  d1$beta <- keep_vec(d1$beta, keep)
  d1$varbeta <- keep_vec(d1$varbeta, keep)
  d1$snp <- keep_vec(d1$snp, keep)
  d1$position <- keep_vec(d1$position, keep)
  d1$MAF <- keep_vec(d1$MAF, keep)
  d2$beta <- keep_vec(d2$beta, ord)
  d2$varbeta <- keep_vec(d2$varbeta, ord)
  d2$snp <- keep_vec(d2$snp, ord)
  d2$position <- keep_vec(d2$position, ord)
  d2$MAF <- keep_vec(d2$MAF, ord)

  d1$snp <- as.character(d1$snp)
  d2$snp <- as.character(d2$snp)
  # Force identical SNP naming/order (defensive)
  d2$snp <- d1$snp
  # Force identical positions as well (defensive; coloc expects matching SNP order)
  d2$position <- d1$position

  if (has_flag("--debug")) {
    cat("[DEBUG] names(d1):", paste(names(d1), collapse=","), "\n")
    cat("[DEBUG] names(d2):", paste(names(d2), collapse=","), "\n")
    cat("[DEBUG] d1 lengths:", length(d1$beta), length(d1$varbeta), length(d1$snp), length(d1$position), "\n")
    cat("[DEBUG] d2 lengths:", length(d2$beta), length(d2$varbeta), length(d2$snp), length(d2$position), "\n")
    cat("[DEBUG] identical snp:", identical(d1$snp, d2$snp), "\n")
    cat("[DEBUG] duplicated snp d1:", sum(duplicated(d1$snp)), " d2:", sum(duplicated(d2$snp)), "\n")
    cat("[DEBUG] first snps:", paste(head(d1$snp, 5), collapse=","), "\n")
  }

  # Provide SNP names as vector names as well (defensive for coloc's internal alignment)
  names(d1$beta) <- d1$snp
  names(d1$varbeta) <- d1$snp
  names(d1$position) <- d1$snp
  if (!is.null(d1$MAF)) names(d1$MAF) <- d1$snp
  names(d2$beta) <- d2$snp
  names(d2$varbeta) <- d2$snp
  names(d2$position) <- d2$snp
  if (!is.null(d2$MAF)) names(d2$MAF) <- d2$snp

  res <- coloc::coloc.abf(d1, d2)
  pp <- res$summary
  list(result=res, summary=pp, n_overlap=nrow(m), dataset_id=dataset_id)
}

write_json <- function(path, obj) {
  jsonlite::write_json(obj, path, pretty=TRUE, auto_unbox=TRUE)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package `jsonlite` not installed. Install with install.packages('jsonlite').")
}

locus_id <- get_arg("--locus-id", default=NULL)
locus_ids_csv <- get_arg("--locus-ids", default="")
candidates_tsv <- get_arg(
  "--candidates-tsv",
  default=file.path(proj_dir, "results", "phases_g", "tier1_coloc_candidates.tsv")
)
tier1_tsv <- get_arg("--tier1-tsv", default=file.path(proj_dir, "results", "phases_g", "tier1_science_supplement.tsv"))
top_n <- as.integer(get_arg("--top-n", default="10"))

load_dataset_map <- function(path) {
  m <- list()
  if (!file.exists(path)) return(m)
  t <- read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="")
  if (!all(c("locus_id", "blood_eqtl_dataset_id") %in% colnames(t))) return(m)
  for (i in seq_len(nrow(t))) {
    m[[t$locus_id[i]]] <- t$blood_eqtl_dataset_id[i]
  }
  m
}

dataset_map <- load_dataset_map(candidates_tsv)

targets <- c()
if (nzchar(locus_ids_csv)) {
  targets <- trimws(strsplit(locus_ids_csv, ",", fixed=TRUE)[[1]])
  targets <- targets[nzchar(targets)]
} else if (!is.null(locus_id)) {
  targets <- c(locus_id)
} else if (file.exists(tier1_tsv)) {
  t1 <- read.table(tier1_tsv, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="")
  if (!("symbol" %in% colnames(t1)) || !("blood_eqtl_top_rsid" %in% colnames(t1))) {
    stop("tier1 TSV does not look like a PHASES-G tier1 table")
  }
  # The stronger-genetics locus_id is symbol + "_" + sentinel rsid in our prep pipeline.
  t1$locus_id <- paste0(t1$symbol, "_", t1$blood_eqtl_top_rsid)
  targets <- head(unique(t1$locus_id), top_n)
  if ("blood_eqtl_dataset_id" %in% colnames(t1)) {
    for (i in seq_len(nrow(t1))) {
      dataset_map[[t1$locus_id[i]]] <- t1$blood_eqtl_dataset_id[i]
    }
  }
} else {
  stop("Provide --locus-id or ensure tier1 TSV exists.")
}

empty_row <- function(lid, ds, status, err) {
  data.frame(
    locus_id = lid,
    dataset_id = if (is.null(ds) || length(ds) == 0 || is.na(ds)) "" else as.character(ds),
    status = status,
    error = gsub("\t", " ", err, fixed=TRUE),
    n_overlap = NA_integer_,
    PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
    stringsAsFactors = FALSE
  )
}

out_rows <- list()
for (lid in targets) {
  cat("[COLOC]", lid, "\n")
  ds <- dataset_map[[lid]]
  out_json <- file.path(out_dir, paste0("coloc_", lid, ".json"))
  row <- tryCatch(
    {
      fit <- harmonize_and_coloc(lid, dataset_id = ds)
      s <- as.list(fit$summary)
      obj <- list(
        locus_id = lid,
        dataset_id = if (is.null(fit$dataset_id) || is.na(fit$dataset_id)) "" else fit$dataset_id,
        status = "ok",
        error = "",
        n_overlap = fit$n_overlap,
        summary = s
      )
      write_json(out_json, obj)
      data.frame(
        locus_id = lid,
        dataset_id = obj$dataset_id,
        status = "ok",
        error = "",
        n_overlap = fit$n_overlap,
        PP.H0 = s[["PP.H0.abf"]],
        PP.H1 = s[["PP.H1.abf"]],
        PP.H2 = s[["PP.H2.abf"]],
        PP.H3 = s[["PP.H3.abf"]],
        PP.H4 = s[["PP.H4.abf"]],
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      empty_row(lid, ds, "error", conditionMessage(e))
    }
  )
  out_rows[[length(out_rows) + 1]] <- row
}

if (length(out_rows) > 0) {
  tab <- do.call(rbind, out_rows)
  sum_path <- file.path(out_dir, "coloc_summary.tsv")
  if (has_flag("--merge-summary") && file.exists(sum_path)) {
    old <- tryCatch(
      read.table(sum_path, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char=""),
      error = function(e) NULL
    )
    if (!is.null(old) && nrow(old) > 0 && "locus_id" %in% colnames(old)) {
      new_ids <- unique(as.character(tab$locus_id))
      old <- old[!old$locus_id %in% new_ids, , drop=FALSE]
      tab <- rbind(old, tab)
      tab <- tab[order(tab$locus_id), , drop=FALSE]
    }
  }
  write.table(tab, sum_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote", sum_path, "\n")
}

cat("Done. Output folder:", out_dir, "\n")

