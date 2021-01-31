# sex_lobster
Investigating the influence of fishing on male and female lobster.

## 1. Sampling

<img align="center" height="300" src="01-sampling/sampling.png"></img>

## 2. Outlier detection

We use `pcadapt`to test for local dapatation in *P. elephas*.
We apply a conservative cut-off by selecting only markers in the top 1% of the P-values distribution.
We discover a total of 833 snps putatively under divergent selection out of 83372snps genotyped on 243 individuals.

## 3. Population structure

### Neutral

Using ADMIXTURE and a Discriminant Analysis of Principal Components (DAPC) available in the [adegenet package](https://www.rdocumentation.org/packages/adegenet/versions/2.0.1), we run our script `population_structure_pal.T`. 
We observe complete panmixia, which is likely due to the life history traits of this species (e.g. long pelagic duration phase). 
This result consider all the putatively neutral genetic markers (25230 SNPs).

Then, we investigate the hypothesis of sex-bias dispersal by calculating the degree of relatedness for male and female, separately.
We detect a subtle difference, males being less relative than female, but it turns out that this difference is not significant (Wilcoxon test, P-value = 0.22).

<img align="center" height="600" src="03-population_structure/relatedness.png"></img>

### Adaptive

We decide to look for population structure using the set of SNPs putatively under divergent selection by performing the script `adaptive_structure.R`.
When we look at these highly differentiated regions of the genome, we highlight that 272 SNPs are linked to gender information, meaning that these markers are probably located to sex chromosome (hereafter called sex-linked markers) and in autosome.

<img align="center" height="600" src="03-population_structure/adaptive_structure.png"></img>

## 4. Sex-linked markers

Genetic diversity, which reflects the potential of a species to cope with environmental change, is significantly different if we include this sex-linked markers than if we remove it, meaning that, without gender information (thanks for collecting it!), we may have biased our genetic diversity estimates and then make wrong recommendation for conservation planning.

## 5. Seascape genomics

We remove sex-linked markers and consider each sex separately (female and male) to perform a distance-based analysis available in the package [vegan](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/capscale) with the script `dbRDA-per-sex.R`. 
We find that males seems more influenced by the seascape than females, with seascape features explaining more genomic variation in males than in females (21% was explained by seascape features versus 12% in female). 

<img align="center" height="600" src="04-dbrda/db-rda.png"></img>

## 6. Analysis of genetic diversity 

We uncover a difference in genetic diversity between males and females, females show higher genetic diversity than male, even when we consider non sex-linked markers.
We hypothesize that fishing may have impacted the population differently, with males being more fished than female, which result in a loss of genetic diversity higher for male than for female.

<img align="center" height="600" src="05-genetic-diversity/het_adaptive_sex.png"></img>
