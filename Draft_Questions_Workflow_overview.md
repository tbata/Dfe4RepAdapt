# RepAdaptDFE: design outline for a comparative DFE analysis

Started Oct 19, 2026 · @Thomas Bataillon @Janek Sendrowski

## What is this ?
A first draft of a RepAdapt comparative DFE study, using [PrimateDFE](https://github.com/Sendrowski/PrimateDFE) ([Sendrowski et al. 2026, Genetics](https://doi.org/10.1093/genetics/iyag252)) as the template. It reviews PrimateDFE's choices and what must change for the [RepAdapt](https://github.com/RepAdapt) data. Aims and hypotheses are in `Hypothesis.md`.

## Review of the PrimateDFE pipeline

PrimateDFE is a single Snakemake workflow built on [fastDFE](https://fastdfe.readthedocs.io/en/latest/) (and soon [Ancestree](https://doi.org/10.64898/2026.09.11.750934)). Per species it runs:

1. **VCF prep**: exons only, normalized, biallelic SNPs.
2. **Outgroups**: minimap2 whole-genome alignment, outgroup bases added to the VCF (to be coordinated with the RepAdapt VCFs).
3. **Annotation**: 0-fold and 4-fold sites (fastDFE `DegeneracyAnnotation`), CpG sites removed.
4. **Polarization**: fastDFE `MaximumLikelihoodAncestralAnnotation` (EST-SFS model, K2, 9 subsampled ingroup haplotypes).
5. **SFS**: fastDFE `Parser` with probabilistic projection and polarization, `TargetSiteCounter` for L.
6. **Ne**: Watterson θ of the 4-fold SFS / 4μ.
7. **DFE**: discrete and gamma–exponential, deleterious-only and full, plus joint fits.
8. **Divergence**: fixed differences to the inferred ancestor, fitted jointly with the SFS.
9. **Comparative**: OLS and PGLS (Brownian, Grafen) on a pruned species tree, quantile epistasis test, πN/πS check.

## What transfers to RepAdapt and what does not:

| Component | Verdict for RepAdapt | Note |
| --- | --- | --- |
| fastDFE annotation, polarization, SFS, inference scripts | Reuse almost as-is | Taxon-agnostic; driven by FASTA, GFF, VCF |
| minimap2 outgroup alignment | Reuse, retune | Plant outgroups are more divergent |
| Divergence counting and joint SFS+divergence fit | Reuse | |
| PGLS regressions | Reuse | With a RepAdapt species tree |

## Data

Inputs per species are the RepAdapt deliverables: realigned BAMs, the reference FASTA and gene GFF, and the VCF from [nextflow\_snp\_calling](https://github.com/RepAdapt/nextflow_snp_calling_linux). That VCF is called with `bcftools mpileup | bcftools call -mv` and filtered only on `AC=AN || MQ<30`, so it is variants-only and minimally filtered — both might matter for DFE work.

**Inclusion criteria :**

- $C1$: ≥ 4 diploid individuals (n = 8 haplotypes) from one genetic cluster; more is better, since RepAdapt samples are often larger than the primate panel.
- $C2$: High-coverage data only. Low-coverage species from the [ANGSD pipeline](https://github.com/RepAdapt/ANGSD_Snakemake_Pipeline) (< 3X) at most enter a folded, deleterious-only side analysis via `realSFS`.
- $C3$: A GFF complete enough to assign degeneracy (e.g. BUSCO completeness threshold on the annotation, to be set).
- $C4$: At least two outgroup assemblies alignable to the reference (see ancestral allele section).

- **Population choice.** One cluster per species from the sNMF output of the [population\_structure](https://github.com/RepAdapt/population_structure) pipeline, since pooling a whole species range mixes demes. 2–3 populations and the pooled range as robustness checks.
- **Genotype filters**, which the RepAdapt VCF lacks: per-genotype DP window (e.g. 0.5× to 2× sample median), GQ ≥ 20, site missingness ≤ 10% in the chosen cluster, biallelic SNPs only.

## Degeneracy annotation

As in PrimateDFE: fastDFE `DegeneracyAnnotation` on each species' own reference and GFF, keeping 0-fold (selected) and 4-fold (neutral) sites. Species with low-quality GFFs are dropped (C3).

Plant annotations first need (CHECK WITH SAM):

- One transcript per gene (longest or canonical isoform), so no site gets two degeneracy classes.
- Matching GFF and FASTA contig names, a frequent failure according to the RepAdapt calling README.
- For stitched references, annotation on the coordinates used in the VCF?

## Ancestral allele annotation and outgroups

Polarize with Ancestree local-tree mode, with fixed-tree mode (EST-SFS model, as in PrimateDFE) as a robustness check. Outgroups are aligned whole-genome with minimap2, keeping one-to-one blocks to avoid paralogs from whole-genome duplications. Species without usable outgroups are dropped (C4), as are those whose divergence estimates break down.

## SFS parsing

fastDFE `Parser` with probabilistic projection and polarization, as in PrimateDFE. The common n is set from the trade-off between species retained and samples per species.

- **Target sites L.** As in PrimateDFE, with the 0-fold/4-fold ratio taken from each species' reference FASTA (`TargetSiteCounter`). L matters mostly for Ne; DFE fractions are invariant to uniform rescaling of both spectra (paper, Appendix Fig. A16).
- **Context filters.** Exclude CpG sites (SNPs and L), as in PrimateDFE. Plants also methylate CHG contexts, so a variant additionally drops CHG sites.
- **gBGC.** GC-conservative (A↔T, C↔G) SFS as the robustness check.

## Estimating Ne

Ne = θ̂W / (4μ) from the 4-fold SFS, as in PrimateDFE. The weak link is μ: plants have few direct estimates and generation times from annual to centuries.

- Compile μ per generation and generation time per species with their source, and flag species where μ is borrowed from a relative.
- Run every regression twice: on log10 Ne and on log10 π4-fold (θ alone, no μ). The second removes μ error from the predictor; the quantile epistasis test still needs μ to unscale, so it inherits that uncertainty.
- For selfing species, θ estimates 4Ne·μ with Ne already reduced by inbreeding; report the selfing rate alongside rather than correcting Ne ad hoc.

## DFE inference

As in PrimateDFE (discrete and gamma–exponential, deleterious-only and full, h = 0.5, ε = 0), so results are comparable across the two studies. In addition:

- The discrete model is the primary report, since the gamma model is expected to saturate at the |S\_d| bound for large-Ne outcrossers such as wind-pollinated trees.
- Joint inference with shared shape within clades such as conifers or Brassicaceae, to test whether close relatives share DFE shape.

## Inference from SFS and divergence counts

Reuse PrimateDFE's joint SFS + divergence model, which turned the ωa–Ne signal significant in primates and is central to question 3.

RepAdapt changes the caveats rather than the method:

- **Independence.** As in the paper, the unit is the reference genome, not the population. RepAdapt maps some species to a relative's reference (17 references for 25 plant species in Whiting et al. 2024), so the effective sample for ωa regressions is the number of distinct references.
- **Saturation.** With distant outgroups, 4-fold divergence saturates and multiple hits bias D\_neut; drop species whose divergence estimates break down (see the polarization section).
- **Time mismatch.** Divergence integrates long-term Ne, whereas the regressor is contemporary Ne; the gap is larger for species with recent bottlenecks or range expansions, which RepAdapt's climate-gradient sampling may capture. (to be discussed with Sam and others)

## Dominance, ploidy and mating system

PrimateDFE's dominance extension (fixed h, co-estimated h, or h decreasing with |S|) can carry over unchanged, but in RepAdapt dominance interacts potentially more crucially because of mating system shifts.

- **Selfing.** In highly selfing species most individuals are homozygous, so their two haplotypes are not independent draws and h has little effect on the SFS. Proposed rule: when F > 0.5, sample one haplotype per individual (n = number of individuals) and treat the h-sensitivity runs as uninformative for that species.

- **Polyploidy.** fastDFE assumes diploid genotypes. Allopolyploids with subgenome-resolved references (to be checked for, e.g., *Trifolium repens*, the ANGSD pipeline's test species) can be analysed per subgenome as diploids; autopolyploids are excluded from the main analysis.

- **Clonality and perenniality** are recorded as covariates, not modelled.

Sensitivity runs: h = 0.5 (main), h = 0.3 and 0.1 fixed, and the recessive h(S) callback already in `infer_dfe.py` (h = 0.4·e^(−0.02|S|)).

## Comparative and regression analysis

OLS and PGLS (Brownian, Grafen) as in PrimateDFE, on a dated RepAdapt species tree.

| Analysis | Response | Predictor(s) | Expected under nearly neutral, no-epistasis model |
| --- | --- | --- | --- |
| Scaled DFE vs Ne | Mass in S bins (−∞,−10\], (−10,−1\], (−1,0\] | log10 Ne (and log10 π4-fold) | Strong bin rises, neutral bin falls |
| Substitution rates | ω, ωna (polymorphism-only) | log10 Ne | Decline with Ne |
| Adaptive rate | ωa, α (SFS + divergence) | log10 Ne | ωa flat or rising; α rising |
| Unscaled DFE invariance | log10 \|s\_q\| at f = 0.4, 0.5, 0.6 | log10 Ne | Slope 0 (equivalently S\_q slope 1) |
| Model-free check | log10 πN/πS | log10 Ne | Slope ≈ −b |
| RepAdapt contrasts | Same responses | log10 Ne + one covariate: selfing, life form, ploidy | Covariate effect = departure from pure Ne rescaling |
| Gene-set DFE (optional) | DFE of repeated-adaptation orthogroups minus background | Within species, then across | No difference if those genes are not unusually constrained |

Report Pagel's λ alongside Grafen's ρ so readers can judge phylogenetic signal.

## Robustness checks

The PrimateDFE checks are kept; the last four rows are new for RepAdapt.

| Check | Contrast | Origin |
| --- | --- | --- |
| Parametrization | Discrete vs gamma–exp, del vs full | PrimateDFE |
| Polarization | Unfolded vs folded; ε fixed vs free | PrimateDFE |
| gBGC | All SNPs vs GC-conservative only (SFS, L and divergence) | PrimateDFE |
| Sample size | n = 8 vs 12 vs 20 | PrimateDFE |
| Dominance | h = 0.5, 0.3, 0.1, h(S) | PrimateDFE |
| Regression model | OLS vs PGLS–Brownian vs PGLS–Grafen | PrimateDFE |
| Predictor | log10 Ne vs log10 π4-fold (no μ) | New |
| Population definition | Single sNMF cluster vs pooled range | New |
| Methylation context | CpG excluded vs CpG + CHG excluded | New |
| GcbGC check? | A,T versus G,C spectra S"W spectra W2S  | New ? or ask Janek |
| Genotype filtering | Default vs stringent DP/GQ windows | New |

Simulation validation (paper Appendix C) should be repeated in SLiM with RepAdapt-relevant settings: partial selfing, range expansion, and n matching the real data.

## Workflow implementation (Janek need yourinput here)

A new Snakemake repository (working name RepAdaptDFE) with PrimateDFE as the template.

- Deliverables: per-species SFS (2–3 populations and the species range) and DFE tables, joinable with the primate data.

```
RepAdaptDFE/
├── resources/                  # species.tsv, samples.tsv, species tree
└── workflow/                   # Snakefile, rules/, scripts/, envs/, profiles/
```


## Open / high leveldecisions and  impacts on analysis etc

- [x] The biggest issue is polarization: if too few RepAdapt species have usable outgroups, questions 2–3 shrink to a folded, deleterious-only analysis.  
  **A:** Species without usable outgroups or with low-quality GFFs are dropped, as in PrimateDFE.
- [x] Should we have a localtree step and have a review of uSFS sensitivity?  
  **A:** Yes. Ancestree local-tree mode is the main polarization, fixed-tree mode the robustness check.
- [x] Should we up the number of individuals?  
  **A:** Set the common n from the trade-off between species retained and samples per species, once per-cluster sample sizes are known.
- [ ] Which RepAdapt species have high-coverage VCFs, ≥ 4 individuals in one cluster, and a usable GFF? (Builds the species table.)
- [x] Do we include the animal species, or keep the first paper plant-only?  
  **A:** Plant-only. Most RepAdapt species are plants (163 of >250), so a plant-only study has the most power and can focus on plant-specific features, which adding animals would dilute.
- [x] Outgroup ladder per species and the 4-fold divergence ceiling for polarization and divergence counts.  
  **A:** No fixed ceiling: drop a species only if its divergence estimates break down or it is a clear regression outlier.
- [x] Source and quality of μ and generation time per species.  
  **A:** π4-fold as the primary predictor; compare Ne regressions across μ from several sources.
- [ ] Selfing and ploidy rules: confirm F > 0.5 threshold and which polyploid references are subgenome-resolved.  
  **A:** F > 0.5 agreed; subgenome resolution still to check.
- [x] Should we point blank exclude clonal species?  
  **A:** No. Clone-correct (one ramet per genet), then include.
- [x] Is the all-sites CDS recall feasible on the RepAdapt BAMs as stored, or do we need the coverage-stats outputs instead?  
  **A:** Not needed: take the 0-fold/4-fold ratio from the reference FASTA, as in PrimateDFE.
- [x] gBGC check: A,T versus G,C spectra (S→W, W→S)?  
  **A:** Stratifying could quantify gBGC strength, but filtering out affected sites (GC-conservative SNPs only) suffices as a robustness check.
- [x] Scope of question 4 (repeated-adaptation orthogroups): this current study or a follow-up?  
  **A:** Tentative run in this study, kept if promising.
- [x] Workflow implementation.  
  **A:** New repository built with PrimateDFE as the template, conda envs per rule, no containers. Species-level inputs in sample sheets; settings varied for contrastive validation as wildcards.


### Sources

- [Sendrowski/PrimateDFE](https://github.com/Sendrowski/PrimateDFE) 
- [Sendrowski et al. 2026, Genetics, iyag252](https://doi.org/10.1093/genetics/iyag252) (for the Methods)
- [RepAdapt GitHub organization](https://github.com/RepAdapt) (SNP calling, ANGSD, population\_structure, highcov\_fst, selective\_sweeps READMEs)
- [Whiting et al. 2024, Nat Ecol Evol](https://pmc.ncbi.nlm.nih.gov/articles/PMC11461274) (as example emblematic of the type of data we have)
