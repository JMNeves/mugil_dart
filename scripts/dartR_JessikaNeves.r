########################
###                  ###
###       DartR      ###
### By Jessika Neves ###
###   November 2020  ###
###                  ###
########################
#
#This script was developed under the dartR version 1.1.6 available at https://cran.r-project.org/src/contrib/Archive/dartR/
#
#Load the package
library(dartR)
#
###Open the dart file###
#
setwd(choose.dir())
#
###The following command opens the csv file and converts it to gl###
dart<- gl.read.dart(filename = "mugil_dart_singlerow.csv", ind.metafile = "metrics_6sp.csv")
##Individual (=specimen or sample) metadata are user specified, and do not come from DArT.
##These metrics are supplied by the user by way of a metafile, provided at the time of inputting the SNP data
##to the genlight object. A metafile is a comma-delimited file, usually named ind_metrics.csv or similar, that
##contains labelled columns. The file must have a column headed id, which contains the individual (=specimen
##or sample labels) and a column headed pop, which contains the populations to which individuals are assigned.
#
#Procducing a heat map with the raw data
dart2<-dart[,order(dart@other$loc.metrics$CallRate, decreasing=TRUE)] #This command order the snps by callrate
gl.plot(dart2, col=c("DeepSkyBlue", "DeepPink1", "Gold"), legend = F)
#
#Plot callrate by snp
gl.report.callrate(dart, plot=T)
gl.report.callrate(dart, method = "ind", plot = TRUE, v = 2)
#
####################################################
##########                                ##########
##########       USEFUL COMMANDS          ##########
##########                                ##########
####################################################
#
nInd(dart)#| returns the number of individuals in the genlight object.
nLoc(dart)#returns the number of loci.
nPop(dart)#returns the number of populations to which the individuals are assigned.
indNames(dart)#returns or sets labels for individuals.
locNames(dart)#returns or sets labels for loci.
alleles(dart)#returns or sets allelic states of each locus for each individual (e.g. ?A/T?).
ploidy(dart)#returns or sets the ploidy of the individuals (normally diploid or 2).
pop(dart)#returns or sets the population to which each individual belongs. Try also levels(pop(gl)) for a list of unique population names.
dart@other$loc.metrics$AvgPIC #The average of the polymorphism information content (PIC) of the Reference and SNP allele rows
dart@other$loc.metrics$PICSnp #The polymorphism information content (PIC) for the SNP allele row
dart@other$loc.metrics$RepAvg # The proportion of technical replicate assay pairs for which the marker score is consistent.
NA.posi(dart)#returns the loci with missing values, that is, loci for which a sequence tag failed to amplify for each individual.
chr(dart)#returns or sets the chromosome for each locus.
position(dart)#returns or sets the position of each SNP in the sequence tag of each locus.
other(dart)#returns of sets miscellaneous information stored as a list.
glSum(dart)#counts the frequency of the alternate allele for each locus.
glNA(dart)#counts the number of missing values for each locus.
glMean(dart)#computes the relative frequency (a proportion in the range 0-1) of the second allele for each locus.
glVar(dart)#computes the variance of the allele frequency distribution for each locus.
gl.report.repavg(dart)#reports reproducibility. SNP datasets generated by DArT have in index, RepAvg, generated by reproducing the data
gl.report.heterozygosity(dart)#Calculates the observed heterozygisities by population from a genlight object and plots as a barchart ordered on heterozygosity.
#independently for 30 RepAvg is the proportion of alleles that give a reproducible result, averaged over both alleles for each locus.
gl.report.pa.pop(dart) #private alleles in each population
#The base frequencies and transition and transversion ratios can in any case be obtained using
gl.report.bases(dart)
glDotProd(dart)#computes the dot products between all pairs of individuals, with centering and scaling.
#
#You can check the names of all available loc.metrics via:
names(dart@other$loc.metrics)
##
####################################################
##########                                ##########
##########     SUBSETTING THE DATA        ##########
##########                                ##########
####################################################
##
#table on individuals per population
table(pop(dart))
#It is easy to create a barplot on the number of individuals per population:
barplot(table(pop(dart)), las=2)
###
####Creating one dataset per popualtion####
#the command mono.rm is very important, as it keeps the invariant (fixed) snps in the array
#Mugil liza
liza <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("liza_SA"))
#checking the data
liza
#Plotting the callrate per snp per individual
gl.report.callrate(liza)
#Mugil curema
curema <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("cure_SA"))
curema
gl.report.callrate(curema)
#Mugil rubrioculus
rubri <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("rubr_SA"))
rubri
gl.report.callrate(rubri)
#Mugil brevirostris
brevi <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("brev_SA"))
brevi
gl.report.callrate(brevi)
#Mugil curvidens from SA
curv_SA <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_SA"))
curv_SA
gl.report.callrate(curv_SA)
#Mugil curvidens from MB
curv_MB <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_MB"))
curv_MB
gl.report.callrate(curv_MB)
#Putting together the two populations of Mugil curvidens
curv_all <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_MB","curv_SA"))
curv_all
gl.report.callrate(curv_all)
#Putting together Mugil curema and M. incilis
cure_inci <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("cure_SA","inci_SA"))
cure_inci
gl.report.callrate(cure_inci)
####
####################################################
##########                                ##########
##########     FILTERING THE DATA         ##########
##########                                ##########
####################################################
##
#A suggested order here is:
#
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
all<-gl.filter.repavg(dart, t=0.97)
all
#Filter to remove all but one of multiple snps in the same fragment
all_sec<-gl.filter.secondaries(all, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
all_sec
gl.report.callrate(all_sec)
#Filter by monomorphic loci - as they do not provide information for population structure and simply slow the analysis)
all_sec_mono<-gl.filter.monomorphs(all_sec, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
all_sec_mono
#
###Other filters###
#
#Calculates the probabilities of agreement with H-W equilibrium based on observed frequencies of reference homozygotes, heterozygotes and
#alternate homozygotes. Uses the exact calculations contained in function prob.hwe() as developed by Wigginton, JE, Cutler, DJ, and Abecasis, GR.
#Arguments:
#alpha:level of significance (per locus) [Default 0.05]
#basis: basis for filtering out loci (any, HWE departure in any one population) [default basis="any"]
#bon:apply bonferroni correction to significance levels for filtering [default TRUE]
#v: verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
all_sec_mono_hwe<-gl.filter.hwe(all_sec_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#
#Attention: the next command takes too long! It can stay overnight!
#Filter out loci with trimmed sequence tags that are too similar (possible paralogues). Only works if
#TrimmedSequence is available in the loci metadata, therefore we use another test data set here.
##Arguments:
#threshold = a threshold Hamming distance for filtering loci [default 0.2]
#rs = number of bases in the restriction enzyme recognition sequence [default = 4]
#pb = switch to output progress bar [default FALSE]
#v = verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
all_sec_mono_hwe_trim<- gl.filter.hamming(all_sec_mono_hwe, threshold=0.2, pb=T, v=5)
#
####Checking the data####
#heat map snps
dart3<-all_sec_mono_hwe_trim[,order(all_sec_mono_hwe_trim@other$loc.metrics$CallRate, decreasing=TRUE)] #ordering the snps by callrate
gl.plot(dart3, col=c("DeepSkyBlue", "DeepPink1", "Gold"), legend = F)
#missing data
gl.report.callrate(all_sec_mono_hwe_trim)
#
##salving the matrix##
gl.write.csv(all_sec_mono_hwe_trim, outfile = "all_sec_mono_hwe_trim.csv")
#
#Additional filters to apply could be for excluding possible loci under selection (gl.outflank)
#checking loci for linkage disequilibrium (gl.report.ld)
#
###########################################################
#####  Filtering by amount of missing data per locus  #####
###########################################################
#
# mono.rm - Remove monomorphic loci [default TRUE]
# recalc - Recalculate the locus metadata statistics if any individuals are deleted in the filtering [default FALSE]
# plot	- specify if a histogram of call rate is to be produced [default FALSE]
#######################
#set 1: 0% missing data
all_0MD <- gl.filter.callrate(all_sec_mono_hwe_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
all_0MD
gl.report.callrate(all_0MD)
#private alleles
gl.report.pa.pop(all_0MD)
#heat map snps
dart4<-all_0MD[,order(all_0MD@other$loc.metrics$CallRate, decreasing=TRUE)] #ordering snps by callrate
gl.plot(dart4, col=c("DeepSkyBlue", "DeepPink1", "Gold"), legend = F)
##salving the matrix
gl.write.csv(all_0MD, outfile = "all_0MD.csv")
#######################
#set 2: up to 20% missing data
all_20MD <- gl.filter.callrate(all_sec_mono_hwe_trim, method = "loc", threshold = 0.8, mono.rm = F, recalc = F, plot = F, v = 2)
all_20MD
gl.report.callrate(all_20MD)
#heat map snps
dart5<-all_20MD[,order(all_20MD@other$loc.metrics$CallRate, decreasing=TRUE)] #ordering snps by callrate
gl.plot(dart5, col=c("DeepSkyBlue", "DeepPink1", "Gold"), legend = F)
##salving the matrix
gl.write.csv(all_20MD, outfile = "all_20MD.csv")
#######################
#set 3: up to 40% missing data
all_40MD <- gl.filter.callrate(all_sec_mono_hwe_trim, method = "loc", threshold = 0.6, mono.rm = FALSE, recalc = F, plot = F, v = 2)
all_40MD
gl.report.callrate(all_40MD)
#heat map snps
dart6<-all_40MD[,order(all_40MD@other$loc.metrics$CallRate, decreasing=TRUE)] #ordering snps by callrate
gl.plot(dart6, col=c("DeepSkyBlue", "DeepPink1", "Gold"), legend = F)
##salving the matrix
gl.write.csv(all_40MD, outfile = "all_40MD.csv")
########################
##It is also possible to filter individuals by amount of missing data
all_40_60 <- gl.filter.callrate(all_40MD, method="ind", threshold = 0.4, mono.rm = FALSE, recalc = FALSE, plot = FALSE, v = 2)
all_40_60
gl.report.callrate(all_40_60)
#
#
####################################################
##########                                ##########
##########           ANALYSES             ##########
##########                                ##########
####################################################
#
################   SET1   ################
#
gl.basic.stats(all_0MD, digits = 4)
##Population structure
#PCoA#
##Check the script for changing the PCoA plot colors
pc_0MD <- gl.pcoa(all_0MD, nfactors=5)
#Plot 1x2
gl.pcoa.plot(pc_0MD,all_0MD, ellipse = FALSE, p = 0.95,labels = "none", hadjust = 1.5, vadjust = 1, xaxis = 1,yaxis = 2)
#Plot 3x4
gl.pcoa.plot(pc_0MD,all_0MD, ellipse = FALSE, p = 0.95,labels = "none", hadjust = 1.5, vadjust = 1, xaxis = 3,yaxis = 4)
#3d
gl.pcoa.plot.3d(pc_0MD,all_0MD)
#The scree plot - The number of dimensions with substantive information content can be determined by examining a scree plot
gl.pcoa.scree(pc_0MD)
##
#structure
gl2structure(all_0MD, outfile = "0MD_struc.str",  exportMarkerNames = F)
#faststructure
gl2faststructure(all_0MD, outfile = "all_0MD_faststructure.str", probar = TRUE)
##
#Neighbor joining
gl.tree.nj(all_0MD, type="phylogram")
#geentic distance
#To calculate Euclidean distances
all_0MD_dist<-gl.dist.pop(all_0MD, method="euclidean")
all_0MD_dist
dendro_0MD<-hclust(all_0MD_dist, method="average")
plot(dendro_0MD)
#
##Genetic variability
#Between = Fst
library(StAMPP)
fst_0MD <-stamppFst(all_0MD, nboots=1, percent=95, nclusters=1)
#nboots-number of bootstraps to perform across loci to generate confidence intervals and p-values
#percent-the percentile to calculate the confidence interval around
#nclusters-the number of proccesor treads or cores to use during calculations.
round(fst_0MD,3)#3=number of digits
#
#Within = heterozygosity
gl.report.heterozygosity(all_0MD)
#
##Species tree
##Output to phylogetic analyses
fasta<-gl2fasta(all_0MD, method=3, outfile="phylo_all_0MD_snp.fasta")
fasta<-gl2fasta(all_0MD, method=4, outfile="phylo_all0MD_semamb_snp.fasta")
#Method 1 - heterozygous positions are replaced by the standard ambiguity codes. The resultant sequence fragments are
#concatenated across loci to generate a single combined sequence to be used in subsequent ML phylogenetic analyses.
#Method=2 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant sequence fragments are concatenated across loci to generate a single composite haplotype to be used in
#subsequent ML phylogenetic analyses.
#Method 3 - heterozygous positions are replaced by the standard ambiguity codes. The resultant SNP bases are concatenated
#across loci to generate a single combined sequence to be used in subsequent MP phylogenetic analyses.
#Method=4 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant SNP bases are concatenated across loci to generate a single composite haplotype to be used in subsequent MP phylogenetic analyses.
#
#SNAPP#
#Convert a genlight object to nexus format suitable for phylogenetic analysis by SNAPP (via BEAUti)
gl2snapp(all_0MD, outfile = "all_0MD_snapp.nex", outpath = tempdir(), v = 2)
#
#####################################################################
#
################   SET2   ################
#
gl.basic.stats(all_20MD, digits = 4)
##Population structure
#PCoA#
pc_20MD <- gl.pcoa(all_20MD, nfactors=5)
#Plot 1x2
gl.pcoa.plot(pc_20MD,all_20MD, ellipse = FALSE, p = 0.95,labels = "pop", hadjust = 1.5, vadjust = 1, xaxis = 1,yaxis = 2)
#Plot 3x4
gl.pcoa.plot(pc_20MD,all_20MD, ellipse = FALSE, p = 0.95,labels = "pop", hadjust = 1.5, vadjust = 1, xaxis = 3,yaxis = 4)
#3d
gl.pcoa.plot.3d(pc_20MD,all_20MD)
#The scree plot - The number of dimensions with substantive information content can be determined by examining a scree plot
gl.pcoa.scree(pc_20MD)
##
#structure
gl2structure(all_20MD, outfile = "20MD_struc.str",  exportMarkerNames = F)
##
#Neighbor joining
gl.tree.nj(all_20MD, type="phylogram")
#genetic distance
#To calculate Euclidean distances
all_20MD_dist<-gl.dist.pop(all_20MD, method="euclidean")
all_20MD_dist
dendro_20MD<-hclust(all_20MD_dist, method="average")
plot(dendro_20MD)
#
##Genetic variability
#Between = Fst
library(StAMPP)
fst_20MD <-stamppFst(all_20MD, nboots=1, percent=95, nclusters=1)
#nboots-number of bootstraps to perform across loci to generate confidence intervals and p-values
#percent-the percentile to calculate the confidence interval around
#nclusters-the number of proccesor treads or cores to use during calculations.
round(fst_20MD,3)#3=number of digits
#
#Within = heterozygosity
gl.report.heterozygosity(all_20MD)
#
##Species tree
##Output for phylogenetic analyses
fasta<-gl2fasta(all_20MD, method=3, outfile="phylo_all_20MD_snp.fasta")
fasta<-gl2fasta(all_20MD, method=4, outfile="phylo_all_20MD_semamb_snp.fasta")
#Method 1 - heterozygous positions are replaced by the standard ambiguity codes. The resultant sequence fragments are
#concatenated across loci to generate a single combined sequence to be used in subsequent ML phylogenetic analyses.
#Method=2 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant sequence fragments are concatenated across loci to generate a single composite haplotype to be used in
#subsequent ML phylogenetic analyses.
#Method 3 - heterozygous positions are replaced by the standard ambiguity codes. The resultant SNP bases are concatenated
#across loci to generate a single combined sequence to be used in subsequent MP phylogenetic analyses.
#Method=4 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant SNP bases are concatenated across loci to generate a single composite haplotype to be used in subsequent MP phylogenetic analyses.
#
#SNAPP#
#Convert a genlight object to nexus format suitable for phylogenetic analysis by SNAPP (via BEAUti)
gl2snapp(all_20MD, outfile = "all_20MD_snapp.nex", outpath = tempdir(), v = 2)
#
#####################################################################
#
################   SET3   ################
#
gl.basic.stats(all_40MD, digits = 4)
##Population structure
#PCoA#
pc_40MD <- gl.pcoa(all_40MD, nfactors=5)
#Plot 1x2
gl.pcoa.plot(pc_40MD,all_40MD, ellipse = FALSE, p = 0.95,labels = "pop", hadjust = 1.5, vadjust = 1, xaxis = 1,yaxis = 2)
#Plot 3x4
gl.pcoa.plot(pc_40MD,all_40MD, ellipse = FALSE, p = 0.95,labels = "pop", hadjust = 1.5, vadjust = 1, xaxis = 3,yaxis = 4)
#3d
gl.pcoa.plot.3d(pc_40MD,all_40MD)
#The scree plot - The number of dimensions with substantive information content can be determined by examining a scree plot
gl.pcoa.scree(pc_40MD)
##
#structure
gl2structure(all_40MD, outfile = "40MD_struc.str",  exportMarkerNames = F)
##
#Neighbor joining
gl.tree.nj(all_40MD, type="phylogram")
#distancia genetica
#To calculate Euclidean distances
all_40MD_dist<-gl.dist.pop(all_40MD, method="euclidean")
all_40MD_dist
dendro_40MD<-hclust(all_40MD_dist, method="average")
plot(dendro_40MD)
#
##Genetic variability
#Between = Fst
library(StAMPP)
fst_40MD <-stamppFst(all_40MD, nboots=1, percent=95, nclusters=1)
#nboots-number of bootstraps to perform across loci to generate confidence intervals and p-values
#percent-the percentile to calculate the confidence interval around
#nclusters-the number of proccesor treads or cores to use during calculations.
round(fst_40MD,3)#3=number of digits
#
#Within = heterozygosity
gl.report.heterozygosity(all_40MD)
#
##Species tree
##Output for phylogenetic analyses
fasta<-gl2fasta(all_40MD, method=3, outfile="phylo_all_40MD_snp.fasta")
fasta<-gl2fasta(all_40MD, method=4, outfile="phylo_all40MD_semamb_snp.fasta")
#Method 1 - heterozygous positions are replaced by the standard ambiguity codes. The resultant sequence fragments are
#concatenated across loci to generate a single combined sequence to be used in subsequent ML phylogenetic analyses.
#Method=2 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant sequence fragments are concatenated across loci to generate a single composite haplotype to be used in
#subsequent ML phylogenetic analyses.
#Method 3 - heterozygous positions are replaced by the standard ambiguity codes. The resultant SNP bases are concatenated
#across loci to generate a single combined sequence to be used in subsequent MP phylogenetic analyses.
#Method=4 - the heterozyous state is resolved by randomly assigning one or the other SNP variant to the individual.
#The resultant SNP bases are concatenated across loci to generate a single composite haplotype to be used in subsequent MP phylogenetic analyses.
#
#SNAPP#
#Convert a genlight object to nexus format suitable for phylogenetic analysis by SNAPP (via BEAUti)
gl2snapp(all_40MD, outfile = "all_40MD_snapp.nex", outpath = tempdir(), v = 2)
#
####################################################
##########                                ##########
##########     ADDITIONAL ANALYSES        ##########
##########                                ##########
####################################################
##
################   Distance Phylogeny on resultant OTUs   ################
##
#The script for distance phylogeny is gl2phylip() which calculates Euclidean distances using dist {stats} then#
#outputs the data in a form suitable for input to the Phylip package written by Joseph Felsenstein
#(http://evolution.genetics.washington.edu/phylip.html) (Felsenstein, 1989). The input file can include replicated
#distance matrices for the purpose of bootstrapping.
phy <- gl2phylip(all_60_40_hwe_trim, outfile="dist_all_60_40_hwe_trim.phy", bstrap=1000)
##
################   SNPRelate   ################
##
#R package SNPRelate is available to undertake principal components analysis and relatedness analysis (Zheng et al., 2012).
gl2gds(all_60_40_hwe_trim, outfile="snprelated.gds")
##
################   AMOVA   ################
##
gl.amova(all_0MD, nperm = 100)
##
################   TREEMIX   ################
##
gl2treemix(dart_5pop, outfile = "treemix_input.gz", outpath = tempdir(),v = 2)
#
####################################################
##########                                ##########
##########     SUMMARY STATISTICS         ##########
##########                                ##########
####################################################
##
#It is possible to calculate summary statistics in arlequin or DNAsp with this dataset
#Make a dataset with fragments instead of snps, without filtering by secondary snps
#Set of data filtered only for reproducibility: all
#Use this set to filter by 0% missing data
all
MD0_inva <- gl.filter.callrate(all, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
MD0_inva
#Then filter by paralogs
MD0_inva_trim<- gl.filter.hamming(MD0_inva, threshold=0.2, pb=T, v=5)
M0D_inva_trim
##saving te matrix##
gl.write.csv(MD0_inva_trim, outfile = "arlequin_0MD_trim.csv")
#Preparing files with the fragments, randomly designating the bases of heterozygous sites
gl2fasta(M0D_inva_trim, method=2, outfile="0MD_inva.fasta")
#It is also possible to create a fasta file with ambiguity codes and use the phase function at DNAsp
gl2fasta(MD0_inva_trim, method=1, outfile="0MD_inva_amb.fasta")
#Use this dataset in arlquin software.
##
####################################################
##########                                ##########
##########  DATASET WITHOUT MUGIL LIZA    ##########
##########                                ##########
####################################################
all_outliza <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("brev_SA","cure_SA", "inci_SA", "rubr_SA", "curv_SA", "curv_MB"))
gl.report.callrate(all_outliza)
#Filter 0% missing data
all_outliza_0MD <- gl.filter.callrate(all_outliza, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
all_outliza_rep<-gl.filter.repavg(all_outliza_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
all_outliza_snp<-gl.filter.secondaries(all_outliza_rep, method="best") 
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis)
all_outliza_mono<-gl.filter.monomorphs(all_outliza_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
all_outliza_hwe<-gl.filter.hwe(all_outliza_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
all_outliza_trim<- gl.filter.hamming(all_outliza_hwe, threshold=0.2, pb=T, v=5)
#
#Structure infile
gl2structure(all_outliza_trim, outfile = "all_outliza_0MD_struc.str",  exportMarkerNames = F)
#
#
#############################################################################################################################
#
####################################################
##########                                ##########
##########     Checking the popualtion    ##########
##########     status after filtering:    ##########
##########                                ##########
####################################################
#
#Mugil liza
liza2 <- gl.keep.pop(all_60_40_hwe_trim, mono.rm = FALSE, pop.list=c("liza_SA"))
liza2
gl.report.callrate(liza2)
#Mugil curema
curema2 <- gl.keep.pop(all_60_40_hwe_trim, mono.rm = FALSE, pop.list=c("cure_SA"))
curema2
gl.report.callrate(curema2)
#Mugil rubrioculus
rubri2 <- gl.keep.pop(all_60_40_hwe_trim, mono.rm = FALSE, pop.list=c("rubr_SA"))
rubri2
gl.report.callrate(rubri2)
#Mugil brevirostris
brevi2 <- gl.keep.pop(all_60_40_hwe_trim, mono.rm = FALSE, pop.list=c("brev_SA"))
brevi2
gl.report.callrate(brevi2)
#Mugil curvidens from SA
curv_SA2 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_SA"))
curv_SA2
gl.report.callrate(curv_SA2)
#Mugil curvidens from MB
curv_MB2 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_MB"))
curv_MB2
gl.report.callrate(curv_MB2)
#Mugil curvindes SA + MB
curv_all2 <- gl.keep.pop(all_0MD, mono.rm = FALSE, pop.list=c("curv_MB","curv_SA"))
curv_all2
gl.report.callrate(curv_all2)
#Mugil curema + M. incilis
cure_inci2 <- gl.keep.pop(all_0MD, mono.rm = FALSE, pop.list=c("cure_SA","inci_SA"))
cure_inci2
gl.report.callrate(cure_inci2)
##
################   Renaming the populations   ################
##
curema<-gl.merge.pop(all_0MD, old=c("cure_SA","inci_SA"), new="curema")
curvidens<-gl.merge.pop(curema, old=c("curv_MB","curv_SA"), new="curvidens")
#creating a file with the new populations name
gl.make.recode.pop(curvidens, outfile = "new_pop_assignments.csv")
gl.make.recode.pop(curema, outfile = "new_pop_assignments2.csv")
#saving it
dart_6pop <- gl.recode.pop(curema, pop.recode="new_pop_assignments2.csv")
dart_5pop <- gl.recode.pop(curvidens, pop.recode="new_pop_assignments.csv")
#checking the data
nPop(dart_6pop) #considering two populations of M. curvidens
nPop(dart_5pop) #considering one popualtion of M. curvidens
nPop(all_0MD)	#considering M. curema and M. incilis separated and the two populations of M. curvidens
##
################   HETROZYGOSITY LEVELS   ################
##
#Observed heterozygosity
gl.report.heterozygosity(dart_5pop)
gl.report.heterozygosity(dart_6pop)
##Expected heterozysosity
#First needs to convert to genind 
dart_5pop_gi<-gl2gi(dart_5pop, v = 1)
dart_6pop_gi<-gl2gi(dart_6pop, v = 1)
#And then calculate with adegenet
library(adegenet)
Hs(dart_5pop_gi)
Hs(dart_6pop_gi)
##
#######################################################################################################
#
#jessika.neves@icbs.ufal.br

