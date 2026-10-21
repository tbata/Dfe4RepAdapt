#!/usr/bin/env Rscript
## rule_chain_dag.R
## Turn Snakemake-style rule chains ("a → b → c") into a workflow DAG figure.
##
## Usage:
##   Rscript rule_chain_dag.R                    # one collapsed "per species" block
##   Rscript rule_chain_dag.R my_dag --expand    # one block per species
##
## Writes <prefix>.dot, plus .svg/.pdf/.png if Graphviz `dot` is on the PATH
## (fallback: DiagrammeR + DiagrammeRsvg [+ rsvg] if installed).

## ---- 1. Inputs: MAnual DAGS

per_species_chain <- "
  filter_genotypes → callable_cds_sites → align_minimap → paf_to_fasta →
  add_outgroups_vcf → annotate_degeneracy → filter_degeneracy → filter_no_cpg →
  annotate_ancestral_allele → compute_sfs → estimate_Ne → infer_dfe →
  count_divergence → infer_dfe_divergence"

cross_species_chain <- "assemble_dfe_table → extract_regression_vars → fit_pgls"

species   <- c("species_A", "species_B", "species_C")  # your species
sfs_sets  <- c("pop1", "pop2", "species_wide")         # SFS flavours per species
new_rules <- c("callable_cds_sites")                   # highlighted as new

## Which rules write / read the deliverables
sfs_writer  <- "compute_sfs"                           # -> results/sfs/
dfe_writers <- c("infer_dfe", "infer_dfe_divergence")  # -> results/tables/dfe/
dfe_reader  <- "assemble_dfe_table"                    # <- results/tables/dfe/

## Optional: real dependencies that the linear chain hides (drawn dotted).
## e.g. count_divergence really needs the annotated alignment, not infer_dfe:

extra_edges <- data.frame(
  from = character(0),   # c("annotate_ancestral_allele")
  to   = character(0)    # c("count_divergence")
)

## ---- 2. Parsing  function

parse_chain <- function(x) {
  x <- gsub("`", "", x, fixed = TRUE)
  x <- gsub("\u2192", "->", x, fixed = TRUE, useBytes = TRUE)  # locale-safe "→"
  nodes <- trimws(unlist(strsplit(x, "->", fixed = TRUE)))
  nodes <- nodes[nzchar(nodes)]
  if (anyDuplicated(nodes)) stop("Duplicated rule in chain: ",
                                 paste(nodes[duplicated(nodes)], collapse = ", "))
  nodes
}

## ---- 3. DOT helpers 
q <- function(x) paste0('"', x, '"')

node_stmt <- function(id, label, fill, shape = "box", style = "rounded,filled") {
  sprintf('    %s [label=%s, shape=%s, style="%s", fillcolor="%s"];',
          q(id), q(label), shape, style, fill)
}

edge_stmt <- function(from, to, attrs = "") {
  sprintf("    %s -> %s%s;", q(from), q(to),
          if (nzchar(attrs)) paste0(" [", attrs, "]") else "")
}

## ---- 4. Build the graph -
build_dot <- function(ps, cs, species, sfs_sets, new_rules, extra_edges,
                      expand = FALSE, rankdir = "TB") {
  stopifnot(sfs_writer %in% ps, all(dfe_writers %in% ps), dfe_reader %in% cs)
  
  col_rule  <- "#EAF0F8"; col_new <- "#FCE3B0"
  col_out   <- "#E3F2E1"; col_cs  <- "#F1E6F5"
  out_edge  <- 'style=dashed, color="#2E7D32"'
  dfe_id    <- "out::dfe"
  units     <- if (expand) species else "per_species"
  
  L <- c("digraph rule_chain {",
         sprintf('  graph [rankdir=%s, fontname="Helvetica", nodesep=0.3, ranksep=0.3];',
                 rankdir),
         '  node  [fontname="Helvetica", fontsize=11];',
         '  edge  [color="#555555", arrowsize=0.7];')
  
  for (u in units) {
    id     <- function(r) paste0(u, "::", r)
    title  <- if (expand) u else sprintf("per species (n = %d)", length(species))
    sp_tag <- if (expand) u else "{species}"
    
    L <- c(L,
           sprintf("  subgraph %s {", q(paste0("cluster_", u))),
           sprintf('    label=%s; style="rounded,dashed"; color="#7A8CA5"; fontsize=12;',
                   q(title)))
    
    for (r in ps) {
      is_new <- r %in% new_rules
      L <- c(L, node_stmt(id(r), if (is_new) paste0(r, "\\n(new)") else r,
                          if (is_new) col_new else col_rule))
    }
    
    sfs_lab <- paste0("results/sfs/", sp_tag, "/\\n",
                      paste(sfs_sets, collapse = " | "))
    L <- c(L,
           node_stmt(id("sfs_out"), sfs_lab, col_out, shape = "folder", style = "filled"),
           edge_stmt(id(head(ps, -1)), id(tail(ps, -1))),
           edge_stmt(id(sfs_writer), id("sfs_out"), out_edge))
    
    if (nrow(extra_edges))
      L <- c(L, edge_stmt(id(extra_edges$from), id(extra_edges$to),
                          'style=dotted, constraint=false'))
    L <- c(L, "  }",
           edge_stmt(id(dfe_writers), dfe_id, out_edge))
  }
  
  ## Shared DFE table directory = hand-off between the two halves
  L <- c(L, node_stmt(dfe_id, "results/tables/dfe/\\n(one table per species)",
                      col_out, shape = "folder", style = "filled"))
  
  cid <- paste0("across::", cs)
  L <- c(L,
         '  subgraph "cluster_across" {',
         '    label="across species"; style="rounded,dashed"; color="#8E6BA8"; fontsize=12;',
         node_stmt(cid, cs, col_cs),
         edge_stmt(head(cid, -1), tail(cid, -1)),
         "  }",
         edge_stmt(dfe_id, paste0("across::", dfe_reader), out_edge),
         "}")
  paste(L, collapse = "\n")
}

## ---- 5. Render functions 

render_dag <- function(dot, prefix) {
  dot_file <- paste0(prefix, ".dot")
  writeLines(dot, dot_file)
  if (nzchar(Sys.which("dot"))) {
    for (fmt in c("svg", "pdf", "png"))
      system2("dot", c(paste0("-T", fmt), if (fmt == "png") "-Gdpi=200",
                       shQuote(dot_file), "-o", shQuote(paste0(prefix, ".", fmt))))
  } else if (requireNamespace("DiagrammeR", quietly = TRUE) &&
             requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
    svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
    writeLines(svg, paste0(prefix, ".svg"))
    if (requireNamespace("rsvg", quietly = TRUE)) {
      rsvg::rsvg_pdf(charToRaw(svg), paste0(prefix, ".pdf"))
      rsvg::rsvg_png(charToRaw(svg), paste0(prefix, ".png"), width = 2000)
    }
  } else {
    message("No renderer found; wrote ", dot_file,
            " (install Graphviz, or DiagrammeR + DiagrammeRsvg).")
  }
  invisible(dot)
}

## ---- 6. Main  
args   <- commandArgs(trailingOnly = TRUE)
expand <- "--expand" %in% args
rest   <- setdiff(args, "--expand")
prefix <- if (length(rest)) rest[1] else "rule_chain_dag"

dot <- build_dot(parse_chain(per_species_chain), parse_chain(cross_species_chain),
                 species, sfs_sets, new_rules, extra_edges, expand = expand)
render_dag(dot, prefix)
message("Done: ", prefix, ".{dot,svg,pdf,png}")

if (interactive() && requireNamespace("DiagrammeR", quietly = TRUE))
  print(DiagrammeR::grViz(dot))

