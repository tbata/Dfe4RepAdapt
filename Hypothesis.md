# Hypotheses

The hypotheses to be tested with the RepAdapt data.

## Aims and questions

RepAdapt spans far deeper divergence than the catarrhines of PrimateDFE: about 300 Myr across the 25 plant species of Whiting et al. 2024 ([Nat Ecol Evol](https://doi.org/10.1038/s41559-024-02514-5)), and the dataset is growing. It also offers a wider range of mating-system and life-history contrasts.

Open points:

* Process all data with the same workflow (probably yes).
* Analyze all data jointly or split by major clade (to be discussed).

### Main questions

1. **Conservation of the scaled DFE.** How much do the masses in S = 4Ne·s bins vary across species, and which features are most stable? This can be tested clade by clade: is the DFE of amino-acid-changing mutations conserved once differences in Ne are accounted for?
2. **Ne rescaling versus landscape change.** Is cross-species variation explained by Ne alone (quantile slope = 1), or does the unscaled DFE shift with Ne, life form, mating system or ploidy?
3. **Adaptive substitution.** Does ωa increase with Ne once divergence counts are added, and does α track Ne only through ωna?
4. **DFE of repeated-adaptation genes.** Do orthogroups repeatedly associated with climate (e.g. the 108 of Whiting et al. 2024) show a different DFE from background genes within the same species? This uses fastDFE's SFS stratification and extends to finer stratifications.

## Tech notes

Dominance becomes a sensitivity analysis, merged with mating system and ploidy (see the design outline). An open question is how to compare species with different mating systems validly, e.g. purely outcrossing against highly selfing.
