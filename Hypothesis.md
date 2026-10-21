We list below a more explicit set of Hypothesis to be tested with the RepAdapt datasets

## Aims and questions

RepAdapt spans far deeper divergence than we have tried on catarrhines (so far about 300 Myr across its 25 plant species, Whiting et al. 2024, [Nat Ecol Evol](https://doi.org/10.1038/s41559-024-02514-5)) and that is fgrowing .. Plus here we can operate with a wider set of  mating-system and life-history contrasts all nested within large kingdom level studies.
One open Qs is how much we want to:   

* process and analyze all data using the samw workflow (probably yes) 
* analyze all data jointly or split in bigger clades (to be discussed)


### Main Qs:

1. **Conservation of the scaled DFE.** How much do the masses in S = 4Ne·s bins vary across species, and which features are most stable?  

The test can be clade by clade emphasizingwhether the DFE of amino-acid-changing mutations is conserved across RepAdapt species once differences in Ne are accounted for.


2. **Ne rescaling versus landscape change.** Is cross-species variation explained by Ne alone (quantile slope = 1), or does the unscaled DFE shift with Ne, life form, mating system or ploidy?  

3. **Adaptive substitution.** Does ωa increase with Ne once divergence counts are added, and does α track Ne only through ωna?  

4. **RepAdapt-specific: DFE of repeated-adaptation genes.** Do  orthogroups repeatedly associated with climate (eg 108 in Whiting et al. 2024) 
show a different DFE from background genes within the same species? 



This uses fastDFE's SFS stratification and can also be extended to finer grained stratifications

## Tech notes

RE mating systems, we propose that the issue of Dominance and DFE becomes a sensitivity analysis here, merged with mating system and ploidy (see that section). 

An open question that remains is how to do valid comparisons for pair of species with different mating system especially from eg purely outcrossing to quite highly selfing?.

