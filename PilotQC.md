# Pilot Genotyping Data QC code

#!/bin/bash

### Input files must only include data for chromosomes 1-22 +/- 23.
### Create directory containing only:  1) the script, 2) the .ped and .map input genotype files, and 3) the 1000 genomes plink binary files after they have undergone QC and had problematic LD regions removed.
### QC parameters employed: sexcheck - keep those >0.8 and <0.2; coverage - keep individuals with >/=98% coverage; call rate - keep SNPs with call rate >/=98%; MAF - keep SNPs with MAF >/=5%; heterozygosity - keep individuals <4sd from mean F; relatedness - keep individuals with IBD </=0.125.


echo "QC pipeline running"

### If data is split per chromosome, combine all chromosomes
cat *.map > GENO1.map
cat *.ped > GENO1.ped

if [ ! -f GENO1.map ]; then echo "Error: No GENO1 created"; fi
if [ -f GENO1.map ]; then echo "GENO1 .ped and .map files created"; fi

if [ ! -f GENO1.ped ]; then echo "Error: No GENO1 created"; fi

### Convert .ped and .map files into plink binary files
/bin/plink2 --file GENO1 --make-bed --out GENO1

if [ ! -f GENO1.bed ]; then echo "Error: No first bed/bim/fam created"; fi
if [ -f GENO1.bed ]; then echo "Plink binary files created"; fi

### Identify individuals with incorrect sex information
plink2 --bfile GENO1 --check-sex --out GENO1sex

if [ ! -f GENO1sex.sexcheck ]; then echo "Error: No sexcheck created"; fi
if [ -f GENO1sex.sexcheck ]; then echo "Sexcheck file created"; fi

### Write individuals with sex check problems to file
/usr/bin/awk '$6<0.8 && $6>0.2 {print $0}' GENO1sex.sexcheck > GENO1_sexerror.txt
if [ -f GENO1_sexerror.txt ]; then echo "Individuals with sex errors identified"; fi
### Remove these individuals
plink2 --bfile GENO1 --remove GENO1_sexerror.txt --make-bed --out GENO2

if [ ! -f GENO2.bed ]; then echo "Error: No GENO2 created"; fi
if [ -f GENO2.bed ]; then echo "GENO2 plink binary files created - sex errors removed"; fi

### Identify and exclude SNPs with poor coverage (<98%) and individuals with poor call rate (<98%)
plink2 --bfile GENO2 --missing --out GENO2_missing

if [ ! -f GENO2_missing.imiss ]; then echo "Error: No GENO2_missing created"; fi
if [ -f GENO2_missing.imiss ]; then echo "Missing SNPs and individuals identified"; fi

### Create plots in R of missingness
Rscript missing.r

if [ ! -f Missingness.png ]; then echo "Error: No first missingness R plots created"; fi
if [ -f Missingness.png ]; then echo "First missingness plots created"; fi


### Write file of SNPs with coverage <98% and remove them from data
/usr/bin/awk '$5>0.02 {print $0}' GENO2_missing.lmiss > GENO2_missingSNPs.txt
plink2 --bfile GENO2 --exclude GENO2_missingSNPs.txt --make-bed --out GENO3

if [ ! -f GENO3.bed ]; then echo "Error: No GENO3 created"; fi
if [ -f GENO3.bed ]; then echo "GENO3 plink binary files created - missing SNPs excluded"; fi

### Re-calculate missing statistics after excluding SNPs with poor coverage and re-plot
plink2 --bfile GENO3 --missing --out GENO3_missing
Rscript missing2.r

if [ ! -f Missingness2.png ]; then echo "Error: No second missingness R plots created"; fi
if [ -f Missingness2.png ]; then echo "Second missingness plots created"; fi

### Write file of individuals with call rate <98% and remove them from data
/usr/bin/awk '$6>0.02 {print $0}' GENO3_missing.imiss > GENO3_missingindividuals.txt
plink2 --bfile GENO3 --remove GENO3_missingindividuals.txt --make-bed --out GENO4

if [ ! -f GENO4.bed ]; then echo "Error: No GENO4 created"; fi
if [ -f GENO4.bed ]; then echo "GENO4 plink binary files created - missing individuals removed"; fi

### Identify and remove SNPs with low MAF
plink2 --bfile GENO4 --freq --out GENO4_freq
usr/bin/awk '$5<0.05 {print $0}' GENO4_freq.frq > GENO4_lowMAFSNPs.txt
plink2 --bfile GENO4 --exclude GENO4_lowMAFSNPs.txt --make-bed --out GENO5 

if [ ! -f GENO4_freq.frq ]; then echo "Error: No GENO4_freq created"; fi
if [ -f GENO4_freq.frq ]; then echo "MAF file created"; fi
if [ ! -f GENO4_lowMAFSNPs.txt ]; then echo "Error: No GENO4_lowMAFSNPs created"; fi
if [ -f GENO4_lowMAFSNPs.txt ]; then echo "Low MAF SNPs identified"; fi
if [ ! -f GENO5.bed ]; then echo "Error: No GENO5 created"; fi
if [ -f GENO5.bed ]; then echo "GENO5 plink binary files created - low MAF SNPs excluded"; fi

### Determine autosomal heterozygosity for each individual
plink2 --bfile GENO5 --het --out GENO5_het

if [ ! -f GENO5_het.het ]; then echo "Error: No GENO5_het created"; fi
if [ -f GENO5_het.het ]; then echo "First heterozygosity file created"; fi

### Plot heterozygosity in R
Rscript heterozygosity.r

if [ ! -f Heterozygosity.png ]; then echo "Error: no first heterozygosity plot created"; fi 
if [ ! -f HET_high.txt ]; then echo "Error: no HET_high file created"; fi
if [ ! -f HET_low.txt ]; then echo "Error: no HET_low file created"; fi
if [ -f Heterozygosity.png ]; then echo "First heterozygosity plots created"; fi

### Remove individuals with outlying heterozygosity
cat HET_high.txt HET_low.txt > HET_outliers.txt
if [ -f HET_outliers.txt ]; then echo "Heterozygosity outliers identified"; fi
if [ ! -f HET_outliers.txt ]; then echo "Error: heterozygosity outliers not identified"; fi
plink2 --bfile GENO5 --remove HET_outliers.txt --make-bed --out GENO6

if [ ! -f GENO6 ]; then echo "Error: no GENO6 files created"; fi
if [ -f GENO6 ]; then echo "GENO6 plink binary files created - heterozygosity outliers removed"; fi

### Re-calculate heterozygosity and plot
plink2 --bfile GENO6 --het --out GENO6_het

if [ ! -f GENO6_het.het ]; then echo "Error: No GENO6_het created"; fi
if [ -f GENO6_het.het ]; then echo "Second heterozygosity file created"; fi

Rscript heterozygosity2.r

if [ ! -f Heterozygosity2.png ]; then echo "Error: no second heterozygosity plot created"; fi
if [ -f Heterozygosity2.png ]; then echo "Second heterozygosity plots created"; fi

### Identify and exclude related pairs of individuals
### LD prune data
plink2 --bfile GENO6 --indep-pairwise 1500 150 0.2 --out GENO6_LDprune

if [ ! -f GENO6_LDprune.prune.in ]; then echo "Error: no prune.in file created"; fi
if [ -f GENO6_LDprune.prune.in ]; then echo "LD pruning - done"; fi

plink2 --bfile GENO6 --extract GENO6_LDprune.prune.in --genome --out GENO6_relatedness

if [ ! -f GENO6_relatedness.genome ]; then echo "Error: no relatedness file created"; fi
if [ -f GENO6_relatedness.genome ]; then echo	"Relatedness file created"; fi 

/usr/bin/awk '$10==1 {print $1,$2,$3,$4,$10}' GENO6_relatedness.genome > GENO6_identical.txt
/usr/bin/awk '{print $3,$4,$5}' GENO6_identical.txt > GENO6_otheridentical.txt

if [ ! -f GENO6_identical.txt ]; then echo "Error: Identical individuals not identified"; fi
if [ -f GENO6_identical.txt ]; then echo "Identical individuals identified"; fi

/usr/bin/awk '$10<1 {print $1,$2,$3,$4,$10}' GENO6_relatedness.genome > GENO6_notidentical.txt

if [ ! -f GENO6_notidentical.txt ]; then echo "Error: non-identical individuals not identified"; fi
if [ -f GENO6_notidentical.txt ]; then echo "Non-identical individuals identified"; fi

/usr/bin/awk '$10>0.125 {print $1,$2,$3,$4,$10}' GENO6_notidentical.txt > GENO6_related.txt

if [ ! -f GENO6_related.txt ]; then echo "Error: Related individuals not identified"; fi
if [ -f GENO6_related.txt ]; then echo "Related individuals identified"; fi

plink2 --bfile GENO6 --remove GENO6_identical.txt --make-bed --out GENO7
plink2 --bfile GENO7 --remove GENO6_otheridentical.txt --make-bed --out GENO8
plink2 --bfile GENO8 --remove GENO6_related.txt --make-bed --out GENO9

if [ ! -f GENO7.bed ]; then echo "Error: no GENO7 plink binary files created"; fi
if [ -f GENO7.bed ]; then echo "GENO7 plink binary files created - identical individuals removed"; fi

if [ ! -f GENO8.bed ]; then echo "Error: no GENO8 plink binary files created"; fi
if [ -f GENO8.bed ]; then echo "GENO8 plink binary files created - other identical individuals removed"; fi

if [ ! -f GENO9.bed ]; then echo "Error: no GENO9 plink binary files created"; fi
if [ -f GENO9.bed ]; then echo "GENO9 plink binary files created - related individuals removed"; fi

### Calculate Hardy-Weinberg Equilibrium statistics for further use if needed
plink2 --bfile GENO9 --hardy --out GENO9_HWE

echo "PIPELINE COMPLETE"

## Rscripts

### missing.r
#!/usr/bin/Rscript
IMISS <- read.table("GENO2_missing.imiss", h=T, as.is=T)
LMISS <- read.table("GENO2_missing.lmiss", h=T, as.is=T)
png("Missingness.png")
oldpar <- par(mfrow = c(1,2))
plot((1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), main ="Ordered individual call rate", xlab="Quantile", ylab="Call rate");
plot((1:dim(LMISS)[1])/(dim(LMISS)[1]-1), sort(1-LMISS$F_MISS), main ="Ordered SNP coverage", xlab="Quantile", ylab="Coverage");
par(oldpar)
dev.off()

### missing2.r
#!/usr/bin/Rscript
IMISS2 <- read.table("GENO3_missing.imiss", h=T, as.is=T)
LMISS2 <- read.table("GENO3_missing.lmiss", h=T, as.is=T)
png("Missingness2.png")
oldpar <- par(mfrow = c(1,2))
plot((1:dim(IMISS2)[1])/(dim(IMISS2)[1]-1), sort(1-IMISS2$F_MISS), main ="Ordered individual call rate", xlab="Quantile", ylab="Call rate");
plot((1:dim(LMISS2)[1])/(dim(LMISS2)[1]-1), sort(1-LMISS2$F_MISS), main ="Ordered SNP coverage", xlab="Quantile", ylab="Coverage");
par(oldpar)
dev.off()

### heterozygosity.r
#!/usr/bin/Rscript

HET <- read.table("GENO5_het.het", h=T, as.is=T)
H <- (HET$N.NM. - HET$O.HOM.) / HET$N.NM.
png("Heterozygosity.png")
oldpar = par(mfrow=c(1,2))
hist(H,50, main="Histogram of Heterozygosity", xlab="Heterozygosity")
abline(v=median(H), lty=3, col='red')
hist(HET$F,50, main="Histogram of F coefficient", xlab="Inbreeding coefficient")
abline(v=median(HET$F), lty=3, col='red')
abline(v=(mean(HET$F) + 4*(sd(HET$F))), lty=3, col='red')
abline(v=(mean(HET$F) - 4*(sd(HET$F))), lty=3, col='red')
par(oldpar)
dev.off()

HEThigh <- HET[HET$F>((mean(HET$F) + 4*(sd(HET$F)))),]
HETlow <- HET[HET$F<((mean(HET$F) - 4*(sd(HET$F)))),]

write.table(HEThigh, file="HET_high.txt", quote=FALSE, row.names=FALSE)
write.table(HETlow, file="HET_low.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

### heterozygosity2.r
#!/usr/bin/Rscript

HET2 <- read.table("GENO6_het.het", h=T, as.is=T)
H2 <- (HET2$N.NM. - HET2$O.HOM.) / HET2$N.NM.
png("Heterozygosity2.png")
oldpar = par(mfrow=c(1,2))
hist(H2,50, main="Histogram of Heterozygosity", xlab="Heterozygosity")
abline(v=median(H2), lty=3, col='red')
hist(HET2$F,50, main="Histogram of F coefficient", xlab="Inbreeding coefficient")
abline(v=median(HET2$F), lty=3, col='red')
par(oldpar)
dev.off()


# Pilot genotyping QC notes

Pilot data:  SNP genotyping data (.ped and .map files) from Illumina Human660-Quad chip – 509 cases, 594398 SNPs.
Software required: Plink.

Step 1:  Convert .ped and .map files to binary files:

./plink --file MS_UKW_ALL_illumina --make-bed --out MS_UKW_ALL_illumina

Problem:  Data comes in two forms (1) combined PED and MAP files for whole genome and (2) separated into data for each chromosome.  On inspection of the MAP files for the combined data it was apparent that the SNPs had incorrect chromosome labeling.  Therefore, using the command line, I combined all MAP files for each chromosome:

Combining all _plink.gen.map files:
cat MS_UKW_0X_illumina_plink.gen.map … MS_UKW_22_illumina_plink.gen.map > ALL_illumina_plink.gen.map

Then used that to create a text file containing 2 columns:

awk ‘{ print $2 $1 }’ ALL_illumina_plink.gen.map > ALL_SNP_Chr_updated.txt

SNP	         Chr
rsUVWXYZ 1

I then used plink (--update-chr) to update/correct the chromosome labeling of the combined data MAP file:

./plink --bfile MS_UKW_ALL_illumina --update-map ALL_SNP_Chr_updated.txt --update-chr --make-bed --out MS_UKW_ALL_illumina_updatedChr

Inspection of the updated MAP file showed that there were still 23 SNPs labeled as Chr 0, 15 labeled as Chr 25 and 136 labeled as Chr 136.

I then checked a random selection of 5 SNPs from each abnormally labeled category (0, 25, 26) on the UCSC genome browser and found:
Chr 0 SNPs:  4/5 mapped to Y chromosome, 1 was not found.
Chr 25 SNPs:  5/5 mapped to both X or Y chromosome.
Chr 26 SNPs:   All  mitochondrial SNPs.

I then used plink to extract data on only the SNPs 1-22 and X(23):

1)	Create txt file of desired SNPs from ALL_SNP_Chr_updated.txt:

awk ‘{ print $1 }’ ALL_SNP_Chr_updated.txt > SNPs_1_23.txt

2) Create bed, bim and fam files of desired SNPs only:

./plink --bfile MS_UKW_ALL_illumina_updatedChr --extract SNPs_1_23.txt --make-bed --out 2011_SNPs_1_23

## Identify individuals with discordant sex information:

The sex chromosomes contain portions that take part in recombination, which are inherited like any autosomal genes, and also portions that are not.  Males are XY therefore cannot be heterozygous for any SNPs that do not lie in the pseudo-autosomal region of the Y chromosome.  Homozygosity for these SNPs in males is expected to be close to 1, and in females is expected to be 0.2.  Therefore, sex discrepancy can be identified by calculating homozygosity rates across all X-chromosome SNPs for each individual and comparing these to expected rates.  Unless individuals with sex discrepancy can be identified and accounted for they should be removed.


./plink –bfile 2011_SNPs_1_23 –check-sex –out 2011_SNPs_1_23

No problems detected:  results written to plink.sexcheck.

## Identify and exclude SNPs with poor coverage:

Missing statistics created for both per individual and per SNP:

./plink --bfile 2011_SNPs_1_23 –missing –out 2011_SNPs_1_23
creates 2 files of missingness:  .imiss (call rate across individuals) and .lmiss (coverage across SNPs).

R used to plot graphs of distribution of call rate per individual and coverage per SNP (see R_plot_GWAS2011_missingness.png).

Then command line used to identify 8196 SNPs with coverage < 98%:

awk ‘$5>0.02 {print $0} 2011_SNPs_1_23.lmiss | wc -l

Then command line used to create txt file of SNPs to be excluded:

awk ‘$5>0.02 {print $0} 2011_SNPs_1_23.lmiss > 2011_missing_SNPs.txt

Plink then used to remove SNPs with poor coverage:

./plink --bfile 2011_SNPs_1_23 --exclude 2011_missing_SNPs.txt --make-bed --out 2011_SNPs_1_23_nomissSNPs

Now missingness statistics were re-calculated after exclusion of SNPs with poor coverage:

./plink --bfile 2011_SNPs_1_23_nomissSNPs --missing --out 2011_SNPs_1_23_nomissSNPs

And new R plot created of call rate per individual and coverage per SNP (see R_plot_GWAS2011_nomissSNPs.png

## Identify and remove individuals with poor call rate:

Command line used to identify 3 individuals with call rate < 98%:

awk ‘$6<0.02 {print $0} 2011_SNPs_1_23_nomissSNPs.imiss | wc -l

Then command line used to create txt file of individuals to be excluded:

awk ‘$6<0.02 {print $0} 2011_SNPs_1_23_nomissSNPs.imiss > 2011_missing_individuals.txt

Plink used to remove 3 individuals with low call rate:

./plink --bfile 2011_SNPs_1_23_nomissSNPs --remove 2011_missing_individuals.txt --make-bed --out 2011_1_23_nomissSNPsIndiv

## Identify and remove SNPs with low MAF:

Plink used to creat a .frq file that contains MAFs:

./plink --bfile 2011_1_23_nomissSNPsIndiv --freq --out 2011_1_23_nomissSNPsIndiv_freq

Command line used to create a txt file of all (35739) SNPs with MAF < 0.01:

awk ‘$5<0.01 {print $2}’ 2011_1_23_nomissSNPsIndiv_freq.frq > 2011_lowMAFSNPs.txt

Plink used to exclude 35739 SNPs with low MAFs:

./plink --bfile 2011_1_23_nomissSNPsIndiv --exclude 2011_lowMAFSNPs.txt --make-bed --out 2011_NomissHighMAF

## Determine autosomal heterozygosity for each individual:

./plink --bfile 2011_NomissHighMAF --het --out 2011_NomissHighMAF_het

Used R to create histograms of heterozygosity and F coefficients (see 2011_heterozygosity.png).

Identified 4 outlier individuals.

Used command line to identify outliers:

awk '$6<-0.04 || $6>0.04 {if (NR!=1) {print $0}}' 2011_NomissHighMAF_het.het > 2011_F_outliers.txt

Used plink to remove these outliers:
./plink --bfile 2011_NomissHighMAF --remove 2011_F_outliers.txt --make-bed --out 2011_NomissHighMAFF

Re-ran heterozygosity measurements:

./plink --bfile 2011_NomissHighMAFF --het --out 2011_NomissHighMAFF_het

Re-ran R code to view new histograms (see 2011_heterozygosity2.png).

## Identify and exclude related pairs of individuals:

First step is to LD prune SNPs (to identify and list only ‘independent’ SNPs - R2<0.2):

./plink --bfile 2011_NomissHighMAFF --indep-pairwise 1500 150 0.2 --out 2011_NomissHighMAFF_thin

Calculates r2 for all SNPs within window of 1500 SNPs then moves the window along 150 SNPs and repeats across whole genome.

Writes ‘independent’ SNPs to .prune.in file and ‘correlated’ SNPs to .prune.out file.

Then run --genome command to calculate genetic relatedness of all individuals:

Identity by state (IBS) is the metric used:  based on average proportion of alleles shared in common at genotyped independent (r2<0.2) SNPs between each pair of individuals in cohort.  Population mean IBS varies depending on the allele frequencies in the population.  The more related two individuals are the greater the proportion of shared alleles.
Identity by descent (IBD) = the degree of recent shared ancestry and can be estimated using genome-wide IBS data:  IBD >0.98 implies duplicate/monozygotic twins, IBD ~0.5 implies first degree relative, IBD ~0.25 implies second degree relative and IBD ~0.125 implies third degree relative – samples excluded if IDB > 0.125.

./plink --bfile 2011_NomissHighMAFF --extract 2011_NomissHighMAFF_thin.prune.in --genome --out 2011_NomissHighMAFF_genome

Then create file of related individuals (PI_HAT score > 0.125 - at least 3rd degree relatives):

awk ‘$10>0.125 {print $1, $2, $3, $4, $10}’ 2011_NomissHighMAFF_genome.genome > 2011_IBD_Related

Found 6 pairs of individuals with high IBD (see 2011_IBD_Related)

Then remove one of each pair of related individuals:

./plink --bfile 2011_NomissHighMAFF --remove 2011_IBD_Related --make-bed --out 2011_NomissHighMAFFunrelated

To remove other member of pair with IBD 1.0 I created a file with just those two individuals (2011_IBD_remove) and used plink to remove them – because might be wrong sample:

./plink --bfile 2011_NomissHighMAFFunrelated --remove 2011_IBD_remove --make-bed --out 2011_NomissHighMAFFunrelated2

## Calculate fit of each SNP with Hardy-Weinberg Equilibrium:

./plink --bfile 2011_NomissHighMAFFunrelated2 --hardy --out 2011_HWE

## Perform principal component analysis to identify ancestral outliers and hidden population structures:

Principal components analysis (PCA) used to identify population (ancestral) outliers.  PCA is multivariate statistical method that produces a number of uncorrelated variables (principal components) that account for the variance in a data matrix containing observations across a number of potentially correlated variables.  The observations are the individuals and the potentially correlated variables are the SNPs.    A principal component model is built using pruned genome-wide genotyping data from populations of known ancestry.  To detect large-scale (continental resolution) ancestral differences one could use the HapMap genotype data from Europe (CEU), Asia (CHB and JPT) and Africa (YRI).  The top 2 principal components are sufficient to separately cluster individuals from the 3 populations and these are applied to the test cohort, which clusters them alongside the HapMap individuals.

Identify SNPs common to 1000 genomes data and my data:

1000 genomes data available on ROCKS at /export/home/SHARED
Copied 1000 genomes binary files (All populations:  g1000_All_Populations_QC_No_Bad_LD and separate populations) and latest version of my data (2011_NomissHighMAFFUnrelated2) to my home directory on ROCKS.


Created text files of SNPs for each dataset:

awk ‘{ print $2 }’ 2011_NomissHighMAFFunrelated2.bim > 2011_SNPs.txt

awk ‘{ print $2 }’ g1000_All_Populations_QC_No_Bad_LD.bim > g1000_SNPs.txt 


Created text file of SNPs common to both 1000 genomes dataset and 2011 dataset:

awk ‘NR==FNR{c[$1]++;next};c[$1] > 0’ g1000_SNPs.txt 2011_SNPs.txt > PCA_SNPs.txt


Used plink to create binary files including only SNPs common to both datasets:

plink --noweb –bfile g1000_All_Populations_QC_No_Bad_LD --extract PCA_SNPs.txt --make-bed --out g1000_PCA

plink --noweb –bfile 2011_NomissHighMAFFunrelated2 --extract PCA_SNPs.txt --make-bed --out 2011_PCA


Use plink to update SNP positions so that they match in 2011 and g1000 data:

awk ‘{print $2, $4}’ g1000_PCA.bim > g1000_SNPpos.txt

plink --noweb --bfile 2011_PCA --update-map g1000_SNPpos.txt --make-bed --out 2011_PCAupdated


Then merge binary files for 2011 and g1000 datasets:

Merger failed due to 155511 mismatching SNPs in terms of alleles (saved to 2011_g1000.missnp) – need to rule out strand flipping:

plink --noweb --bfile 2011_PCAupdated --flip 2011_g1000.missnp --make-bed --out 2011_PCAupdatedflip

plink --noweb --bfile 2011_PCAupdatedflip --bmerge g1000_PCA.bed g1000_PCA.bim g1000_PCA.fam --make-bed --out 2011flip_g1000


Merger failed, but this corrected the mismatching for all except 276 SNPs, which need flipping back to original configuration for inspection/removal:

plink --noweb --bfile 2011_PCAupdatedflip --flip 2011flip_g1000.missnp --make-bed --out 2011_PCAupdatedflipcorr


Tried merger again:

plink --noweb --bfile 2011_PCAupdatedflipcorr --bmerge g1000_PCA.bed g1000_PCA.bim g1000_PCA.fam --make-bed --out 2011flipcorr_g1000

Failed as 276 SNPs still mismatching in terms of alleles.


Used plink to extract binary files of mismatching SNPs only from 2011 and g1000 datsets:

plink --noweb --bfile 2011_PCAupdatedflipcorr --extract 2011flip_g1000.missnp --make-bed --out 2011_PCAflipmismatches

plink --noweb --bfile g1000_PCA --extract 2011flip_g1000.missnp --make-bed --out g1000_PCAflipmismatches


Visual inspection of bim files showed that all mismatched SNPs have one allele assigned N in 2011 data, but two alleles in g1000 data.

Removed mismatching SNPs from both 2011 and g1000 datasets:

plink --noweb --bfile 2011_PCAupdatedflipcorr --exclude 2011flip_g1000.missnp --make-bed --out 2011_PCAupdatedflipcorrex

plink --noweb --bfile g1000_PCA --exclude 2011flip_g1000.missnp --make-bed --out g1000_PCAex


Merged datasets:

plink --noweb --bfile 2011_PCAupdatedflipcorrex --bmerge g1000_PCAex.bed g1000_PCAex.bim g1000_PCAex.fam --make-bed --out 2011_g1000_PCA


Removed LD regions known to disrupt PCA:

awk '($1==1)&&($4>=48000000)&&($4<=52000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==2)&&($4>=86000000)&&($4<=100500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==2)&&($4>=134500000)&&($4<=138000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==2)&&($4>=183000000)&&($4<=190000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==3)&&($4>=47500000)&&($4<=50000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==3)&&($4>=83500000)&&($4<=87000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==3)&&($4>=89000000)&&($4<=97500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==5)&&($4>=44500000)&&($4<=50500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==5)&&($4>=98000000)&&($4<=100500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==5)&&($4>=129000000)&&($4<=132000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==5)&&($4>=135500000)&&($4<=138500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==6)&&($4>=25500000)&&($4<=33500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==6)&&($4>=57000000)&&($4<=64000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==6)&&($4>=140000000)&&($4<=142500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==7)&&($4>=55000000)&&($4<=66000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==8)&&($4>=8000000)&&($4<=12000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==8)&&($4>=43000000)&&($4<=50000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==8)&&($4>=112000000)&&($4<=115000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==10)&&($4>=37000000)&&($4<=43000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==11)&&($4>=46000000)&&($4<=57000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==11)&&($4>=87500000)&&($4<=90500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==12)&&($4>=33000000)&&($4<=40000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==12)&&($4>=109500000)&&($4<=112000000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt
awk '($1==20)&&($4>=32000000)&&($4<=34500000){print $2}'  2011_g1000_PCA.bim >> badLDregionSNPs.txt 


LD pruned merged dataset using ROCKS:

nano LDprune.sh:

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -b n
#$ -wd /home/wpmjh18/PCA/input/prune
#$ -l h_vmem=10G
#$ -l h_rt=30:00:00

plink --noweb --bfile 2011_g1000_PCA --indep-pairwise 1500 150 0.2 --out 2011_g1000_PCA


chmod 744 james_LDprune

qsub LDprune.sh


Then created bed, bim and fam files including only pruned SNPs for combined dataset:

plink --noweb  --bfile 2011_g1000_PCA --extract 2011_g1000_PCA.prune.in --make-bed --out 2011_g1000_PCApruned

Then create .par file to run EIGENSTRAT:

vi 2011_g1000_PCA.par

genotypename: 2011_g1000_PCApruned.bed
snpname: 2011_g1000_PCApruned.bim
indivname: 2011_g1000_PCApruned.fam
evecoutname: 2011_g1000_PCApruned.pca.evec
evaloutname: 2011_g1000_PCApruned.eval
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 0
qtmode: 0

Submitted the job to ROCKS:

vi eigenstrat.sh
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q all.q
/share/apps/EIG-6.1.4/bin/smartpca -p 2011_g1000_PCA.par
qsub eigenstrat.sh

The .pca.evec output files contains principal components that were plotted in R

Excluded MS individuals that lie 2 standard deviations away from the mean of first two principal components of combined MS and g1000 data:

Read combined g1000 and MS dataset (2011_g1000_PCApruned.pca.evec) into R.
Label MS patients and g1000 individuals within combined file as previous (using fam files, binding and merging).
Write out a list of only MS individuals from combined dataset (combPCA) - using subset.
Write out a list of all cases who exceed 2SDs from the mean of PC1 of MS cases only.
Repeat to create 4 files in total:  PC1toohigh, PC1toolow, PC2toohigh, PC2toolow.  NB:  Ensure that the mean and SD used is that of the original MS only dataset (don’t re-calculate mean and SD after each exclusion).
Exclude PC outliers as described in file PC5 to create a ‘cutMS’ dataset.


# Second attempt of PCA using updated genome build

First repeat attempt did not extract only SNPs common to both datasets first - strange PCA plot produced.
Therefore, repeated the analysis after selecting only SNPs common to both datasets before merger.  First repeat attempt flagged 42 SNPs that required strand flipping - flipped these before proceding with second repeat attempt.

plink --noweb --bfile 2011_NomissHighMAFFunrelated2-updated-AllChrom --flip 2011newbuild_g1000.missnp --make-bed --out 2011_NomissHighMAFFunrelated2-updated-AllChrom_flip


Created text files of SNPs for each dataset:

awk ‘{print $2}’ 2011_NomissHighMAFFunrelated2-updated-AllChrom_flip.bim > 2011newbuild_SNPs.txt

awk ‘{print $2}’ g1000_All_Populations_QC_No_Bad_LD.bim > g1000_SNPs.txt 


Created text file of SNPs common to both g1000 dataset and MS dataset:

awk ‘NR==FNR{c[$1]++;next};c[$1] > 0’ g1000_SNPs.txt 2011newbuild_SNPs.txt > MSg1000newbuild_SNPs.txt


Used plink to create binary files including only SNPs common to both datasets:

plink --noweb --bfile g1000_All_Populations_QC_No_Bad_LD --extract MSg1000newbuild_SNPs.txt --make-bed --out g1000newbuild

plink --noweb --bfile 2011_NomissHighMAFFunrelated2-updated-AllChrom_flip --extract MSg1000newbuild_SNPs.txt --make-bed --out 2011newbuild


Merged MS data with g1000 data:

plink --noweb --bfile 2011newbuild --bmerge g1000newbuild.bed g1000newbuild.bim g1000newbuild.fam --make-bed --out 2011newbuild_g1000

Succeeded, but with warnings of different SNP positions for 77 SNPs.


Use plink to update SNP positions so that they match in 2011 and g1000 data:

awk ‘{print $2, $4}’ g1000newbuild.bim > g1000_SNPpos.txt

plink --noweb --bfile 2011newbuild --update-map g1000_SNPpos.txt --make-bed --out 2011newbuild_updatedPos


Repeat merger:

plink --noweb --bfile 2011newbuild_updatedPos --bmerge g1000newbuild.bed g1000newbuild.bim g1000newbuild.fam --make-bed --out 2011newbuild_g1000


Then LD pruned combined dataset:

nano LDprune.sh:

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -b n
#$ -wd /home/wpmjh18/PCA/input/prune
#$ -l h_vmem=10G
#$ -l h_rt=30:00:00

plink --noweb --bfile 2011newbuild_g1000 --indep-pairwise 1500 150 0.2 --out 2011newbuild_g1000


chmod 744 james_LDprune

qsub LDprune.sh

Created .prune.in (104504 SNPs) and .prune.out (204500 SNPs) files.

Then created bed, bim and fam files including only pruned SNPs for combined dataset:

plink --noweb  --bfile 2011newbuild_g1000 --extract 2011newbuild_g1000.prune.in --make-bed --out 2011newbuild_g1000pruned


Then create .par file to run EIGENSTRAT:

nano 2011newbuild_g1000.par

genotypename: 2011newbuild_g1000pruned.bed
snpname: 2011newbuild_g1000pruned.bim
indivname: 2011newbuild_g1000pruned.fam
evecoutname: 2011newbuild_g1000pruned.pca.evec
evaloutname: 2011newbuild_g1000pruned.eval
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 0
qtmode: 0


Wrote script to submit the job to ROCKS:

nano eigenstrat.sh

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -b n
#$ -wd /home/wpmjh18/PCA/input/prune
#$ -q all.q

/share/apps/EIG-6.1.4/bin/smartpca -p /home/wpmjh18/PCA/bin/2011newbuild_g1000.par


Submitted job to ROCKS:

qsub eigenstrat.sh


The .pca.evec output file contains principal components that were plotted in R:
1) PC1 v PC2 of MS and g1000 data combined.
2) PC1 v PC2 of MS data extracted from combined dataset.
3) MS patients lying outside of 2sds of the mean of PC1 and PC2 were excluded.
4) Final MS patient cohort plotted alone.



# Final numbers:

509 individuals initially:
	3 individuals removed due to call rate <98%
	4 individuals removed due to outlying heterozygosity
	4 individuals removed due to sample duplication
	4 individuals removed due to cryptic relatedness
	5 individuals removed due to outlying ancestry (PCA)
489 individuals included

594398 SNPs
	174 SNPs excluded as unable to update chromosome label
8196 SNPs excluded due to coverage <98%
	35739 SNPs excluded due to MAF <0.01
550289 SNPs included

There are 2 individuals I excluded that the GWAS included because they are potential duplicates (contained in the 4 patients removed due to duplication - didn’t matter to GWAS as both confirmed cases, but matters to me as the clinical data must match the genotyping data).

There are 10 more individuals excluded by the IMSGC GWAS that I have included:
-	3 due to SNP QC
-	1 due to gender
-	6 unexplained

SNP QC cases:
89807_A11:
Missing:  8269 SNPs /  594224 SNPs:   0.01392 missingness (call rate 98.618%).

89800_D02:
Missing:    6330 SNPs / 594224 SNPs:  0.01065 missingness (call rate 98.935%).

89799_H08:  
Missing:  11303 SNPs / 594224 SNPs:  F_miss = 0.01902  - We used a cut off of 98%, which included this individual with a call rate of 98.108%.  GWAS used 99% cut off.

Gender case:
95328_B08: Plink sex-check F value = 1 = definitely male (X heterozygosity), so unclear why this patient was excluded from GWAS.

Unexplained cases – not present in files of excluded individuals, but says “Out” on spreadsheet:
89806_A04
116278_H03
116731_D09
92028_A03
95328_H04
89807_A11

None of these samples appear in either of the GWAS files listing the excluded individuals.  Supplementary notes of the GWAS says that they analysed 481 of our cases, so these cases were probably removed.  The GWAS used different methods to measure heterozygosity and gender, but unlikely to have been due to that.  The GWAS file for those excluded during PCA lists only 2 individuals, which is surprising.  Also, the GWAS used HapMap phase 2 individuals for the PCA and included cases from multiple genetic backgrounds, whereas I have used the more recent 1000 genomes data and the cohort is relatively genetically homogenous–this may explain the differences.

# Imputation
For imputation the data must be of the same genome build as the chosen reference.  For my data build 37 is required (aligned to HRC dataset).  For the PCA (above) I used plink to update the SNP positions and strands, but there were 276 SNPs that mismatched due to different alleles (N) when compared with the reference genome.  In order to correct for these and ensure most accurate build update use the perl script recommended by the Michigan Imputation website:

http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip

Using the updated SNP information contained in this file:

ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

And using a file containing MAFs for each SNP in dataset created using plink:

plink --noweb --bfile 2011_NomissHighMAFFunrelated2 --freq --out 2011_NomissHighMAFFunrelated2.frq

The perl script outputs file summarising the differences between my data and the reference HRC dataset in terms of SNP positions, alleles and strand.  The script also outputs a script of cutomised plink commands to run to update my data (Run-plink.sh) – change it to executable, add '--noweb' to start of each command and run the sript:

chmod -x Run-plink.sh
./Run-plink.sh

This script outputs updated dataset in whole form and by chromosome.

Further prepare data for upload:

Change chromosome labelling to format recognised by Michigan imputation server:

First convert bed, bim and fam files to ped and map files:

plink2 --bfile 2011_NomissHighMAFFunrelated2-updated-chr1 --recode --out 2011_NomissHighMAFFunrelated2-updated-Chrom1

Then change chromosome labelling and convert to vcf file:

plink2 --ped 2011_NomissHighMAFFunrelated2-updated-Chrom1.ped --map 2011_NomissHighMAFFunrelated2-updated-Chrom1.map --output-chr M --recode vcf-iid --out 2011_NomissHighMAFFunrelated2-updated-Chrom1

Now sort and bgzip file:

./vcf-sort 2011_NomissHighMAFFunrelated2-updated-Chrom1.vcf | bgzip -c > 2011_NomissHighMAFFunrelated2-updated-Chrom1.vcf.gz

Repeat for all chromosomes.

Now upload data to Michigan imputation server.
N.B. Chose Eagle for phasing of autosomal chromosomes, but did not allow this for X chromosome, so chose shapeIT for the X chromosome.
X chromosome imputation failed as imputation pipeline is not yet functional.


# QC of imputed data

1. Convert .vcf files to plink binary files:

nano convert_vcf_plink.sh

for i in {1..22}

do

plink2 --vcf chr${i}.dose.vcf --make-bed --out ch${i}

done

chmod 744 convert_vcf_plink.sh
bash convert_vcf_plink.sh


2. Identify SNPs with R2 (INFO) score >/=0.8:

nano find_

cd ~/genomeQC/pilot/imputed/input

for i in {1..22}

do

awk '$7>0.8 {print $1}' ch${i}.info > chr${i}_SNPs_R2good.txt

done

chmod 744 find_


3. Keep only data for these SNPs:

nano extract_SNPs_R2good.sh

cd ~/genomeQC/pilot/imputed/input

for i in {1..22}

do

plink --noweb --bfile ch${i} --extract chr${i}_SNPs_R2good.txt --make-bed --out chr${i}goodR2

done

chmod 744 extract_SNPs_R2good.sh
bash extract_SNPs_R2good.sh


4.  Exclude SNPs with MAF < 1%:

nano keep_SNPs_goodMAF.sh

cd ~/genomeQC/pilot/imputed/input

for x in {1..22}

do

plink --noweb --bfile chr${x}goodR2 --freq --out chr${x}goodR2

awk ‘$5<0.01 {print $2}’ chr${x}goodR2.frq > chr${x}_lowMAF.txt

plink --noweb --bfile chr${x}goodR2 --exclude chr${x}_lowMAF.txt --make-bed --out chr${x}goodR2MAF

done

chmod 744 extract_SNPs_goodMAF.sh
bash extract_SNPs_goodMAF.sh


5. Repeat relatedness exclusions:

First merge binary files:

plink2 --bfile chr1goodR2MAF --bmerge allchrfiles.txt --make-bed --out allchrgoodR2MAF

where allchrfiles.txt is:
chr2goodR2MAF.bed chr2goodR2MAF.bim chr2goodR2MAF.fam
	...		...		...
	
Then run:

plink2 --bfile allchrgoodR2MAF --genome

Write file of individuals with PIHAT score > 0.125 and view file:

awk ‘$10>0.125 {print $1, $2, $3, $4, $10}’ plink.genome > related.txt
less related.txt

There were no related individuals.


6. Repeat PCA:

(i) Exclude individuals 2sds away from the mean of PC1 and PC2 from PCA when combined with g1000 data:


