# sex_lobster
Investigating the influence of fishing on male and female lobster.

## 1. Sampling

<img align="center" height="300" src="01-sampling/sampling.png"></img>

## 2. Outlier detection

We use `pcadapt`to test for local dapatation in *P. elephas*.
We applied a conservative cut-off by selecting only markers in the top 1% of the P-values distribution.
We discovered a total of 833 snps putatively under divergent selection out of 83372snps genotyped on 243 individuals.

## 3. Population structure

### Neutral

Using ADMIXTURE and a Discriminant Analysis of Principal Components (DAPC) available in the [adegenet package](https://www.rdocumentation.org/packages/adegenet/versions/2.0.1), we run our script `population_structure_pal.T`. 
We found that all indivdiuals sampled are part of the same population, this is probably due to the long pelagic duration phase (up to 4 months). 
This was found considering all the neutral genetic markers we had (25230 SNPs).

Then, we wanted to investigate the hypothesis of sex-bias dispersal by testing the difference in relatedness for male versus female.
We detected a subtle difference but it turns out that it was not significant (using Wilcoxon test).

<img align="center" height="600" src="03-population-structure/relatedness.png"></img>

### Adaptive

We decided to look for population structure using the set of SNPs putatively under divergent selection by performing the script , we run our script `adaptive_structure.T`.
When we look at these highly differentiated region of the genomes, we found that some regions are linked to gender information, meaning that we have found some sex-linked markers.

<img align="center" height="600" src="adaptive_structure.png"></img>

## 4. Sex-linked markers

<img align="center" height="600" src="adaptive_structure.png"></img>

Genetic diversity, which reflects the potential of a species to cope with environmental change, was different if we include this sex-linked markers than if we remove it, meaning that, without gender information (thanks for collecting it!), we may have biased our genetic diversity estimates and then make wrong recommendation for conservation planning.

## 5. Seascape genomics

When we removed these sex-linked markers and considering each sex separately (female and male), we performed a distance-based analysis available in the package [vegan](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/capscale) with the script `dbRDA-per-sex.R`. 
We could find that difference exist between males and females with seascape features explaining more genomic variation in males than in females (21% was explained by seascape features versus 12%). 

<img align="center" height="600" src="04-db-rda/db-rda.png"></img>

## 5. Analysis of genetic diversity 

We uncover a difference in genetic diversity between genders, female showed higher genetic diversity than males, even when we consider markers not correlated with sex.
We hypothesize that fishing may have impacted the population differently, with males being more fished than female, which result in a loss of genetic diversity higher for male than for female.

<img align="center" height="600" src="05-genetic-diversity-sex.png"></img>
