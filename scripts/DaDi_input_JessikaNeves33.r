########################
###                  ###
###   DartR - DADI   ###
### By Jessika Neves ###
###   March 2021     ###
###                  ###
########################
#
#This script was developed under the dartR version 1.1.6 available at https://cran.r-project.org/src/contrib/Archive/dartR/
#
##Producing a dataset to run DaDi
##As the analysis requeires one dataset per species, it is needed to filter by species before the other filters
#
#Load the packages
library(dartR)
library(radiator)
#
###Open the dart file###
#
setwd(choose.dir())
#
###The following command opens the csv file and converts it to gl###
dart<- gl.read.dart(filename = "mugil_dart_singlerow.csv", ind.metafile = "metrics_6sp.csv")
##
############################################
############ Mugil brevirostris ############
############################################
##
##Keeping the population
brev33 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("brev_SA"))
##
##Filtering
##
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
brev33_rep<-gl.filter.repavg(brev33, t=0.97)
brev33_rep
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis
brev33_mono<-gl.filter.monomorphs(brev33_rep, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
brev33_mono
#Filter by paralogs
brev33_trim<- gl.filter.hamming(brev33_mono, threshold=0.2, pb=T, v=5)
brev33_trim
##################
#Filter 0% missing data
brev33_0MD <- gl.filter.callrate(brev33_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
brev33_0MD
#Filter to remove all but one of multiple snps in the same fragment
brev33_0MD_1SNP<-gl.filter.secondaries(brev33_0MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
brev33_0MD_1SNP
##
##Convert to VCF using radiator package
##
library(radiator)
vcf_brevirostris33_0MD<-genomic_converter(brev33_0MD_1SNP, strata = NULL, output = "vcf",
                                        filename = "brevirostris_0MD33", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
##################
#Filter 40% missing data
brev33_40MD <-gl.filter.callrate(brev33_trim, method = "loc", threshold = 0.6, mono.rm = F, recalc = F, plot = F, v = 2)
brev33_40MD
#Filter to remove all but one of multiple snps in the same fragment
brev33_40MD_1SNP<-gl.filter.secondaries(brev33_40MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
brev33_40MD_1SNP
##Convert to VCF using radiator package
##
vcf_brevirostris33_40MD<-genomic_converter(brev33_40MD_1SNP, strata = NULL, output = "vcf",
                                           filename = "brevirostris_40MD33", parallel.core = parallel::detectCores() - 1,
                                           verbose = TRUE)
##
##
############################################
############    Mugil curema    ############
############################################
##
##Keeping the population
cure33 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("cure_SA", "inci_SA"))
##
##Filtering
##
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
cure33_rep<-gl.filter.repavg(cure33, t=0.97)
cure33_rep
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis
cure33_mono<-gl.filter.monomorphs(cure33_rep, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
cure33_mono
#Filter by paralogs
cure33_trim<- gl.filter.hamming(cure33_mono, threshold=0.2, pb=T, v=5)
cure33_trim
##################
#Filter 0% missing data
cure33_0MD <- gl.filter.callrate(cure33_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
cure33_0MD
#Filter to remove all but one of multiple snps in the same fragment
cure33_0MD_1SNP<-gl.filter.secondaries(cure33_0MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
cure33_0MD_1SNP
##
##Convert to VCF using radiator package
##
library(radiator)
vcf_curema33_0MD<-genomic_converter(cure33_0MD_1SNP, strata = NULL, output = "vcf",
                                        filename = "curema_0MD33", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
##################
#Filter 40% missing data
cure33_40MD <-gl.filter.callrate(cure33_trim, method = "loc", threshold = 0.6, mono.rm = F, recalc = F, plot = F, v = 2)
cure33_40MD
#Filter to remove all but one of multiple snps in the same fragment
cure33_40MD_1SNP<-gl.filter.secondaries(cure33_40MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
cure33_40MD_1SNP
##Convert to VCF using radiator package
##
vcf_curema33_40MD<-genomic_converter(cure33_40MD_1SNP, strata = NULL, output = "vcf",
                                           filename = "curema_40MD33", parallel.core = parallel::detectCores() - 1,
                                           verbose = TRUE)
##
##
############################################
############   Mugil curvidens  ############
############################################
##
##Keeping the population
curv33 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_SA", "curv_MB"))
##
##Filtering
##
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
curv33_rep<-gl.filter.repavg(curv33, t=0.97)
curv33_rep
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis
curv33_mono<-gl.filter.monomorphs(curv33_rep, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
curv33_mono
#Filter by paralogs
curv33_trim<- gl.filter.hamming(curv33_mono, threshold=0.2, pb=T, v=5)
curv33_trim
##################
#Filter 0% missing data
curv33_0MD <- gl.filter.callrate(curv33_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
curv33_0MD
#Filter to remove all but one of multiple snps in the same fragment
curv33_0MD_1SNP<-gl.filter.secondaries(curv33_0MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
curv33_0MD_1SNP
##
##Convert to VCF using radiator package
##
library(radiator)
vcf_curvidens33_0MD<-genomic_converter(curv33_0MD_1SNP, strata = NULL, output = "vcf",
                                        filename = "curvidens_0MD33", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
##################
#Filter 40% missing data
curv33_40MD <-gl.filter.callrate(curv33_trim, method = "loc", threshold = 0.6, mono.rm = F, recalc = F, plot = F, v = 2)
curv33_40MD
#Filter to remove all but one of multiple snps in the same fragment
curv33_40MD_1SNP<-gl.filter.secondaries(curv33_40MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
curv33_40MD_1SNP
##Convert to VCF using radiator package
##
vcf_curvidens33_40MD<-genomic_converter(curv33_40MD_1SNP, strata = NULL, output = "vcf",
                                           filename = "curvidens_40MD33", parallel.core = parallel::detectCores() - 1,
                                           verbose = TRUE)
##
##
############################################
############  Mugil rubrioculus ############
############################################
##
##Keeping the population
rubr33 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("rubr_SA"))
##
##Filtering
##
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
rubr33_rep<-gl.filter.repavg(rubr33, t=0.97)
rubr33_rep
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis
rubr33_mono<-gl.filter.monomorphs(rubr33_rep, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
rubr33_mono
#Filter by paralogs
rubr33_trim<- gl.filter.hamming(rubr33_mono, threshold=0.2, pb=T, v=5)
rubr33_trim
##################
#Filter 0% missing data
rubr33_0MD <- gl.filter.callrate(rubr33_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
rubr33_0MD
#Filter to remove all but one of multiple snps in the same fragment
rubr33_0MD_1SNP<-gl.filter.secondaries(rubr33_0MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
rubr33_0MD_1SNP
##
##Convert to VCF using radiator package
##
library(radiator)
vcf_rubrioculus33_0MD<-genomic_converter(rubr33_0MD_1SNP, strata = NULL, output = "vcf",
                                        filename = "rubrioculus_0MD33", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
##################
#Filter 40% missing data
rubr33_40MD <-gl.filter.callrate(rubr33_trim, method = "loc", threshold = 0.6, mono.rm = F, recalc = F, plot = F, v = 2)
rubr33_40MD
#Filter to remove all but one of multiple snps in the same fragment
rubr33_40MD_1SNP<-gl.filter.secondaries(rubr33_40MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
rubr33_40MD_1SNP
##Convert to VCF using radiator package
##
vcf_rubrioculus33_40MD<-genomic_converter(rubr33_40MD_1SNP, strata = NULL, output = "vcf",
                                           filename = "rubrioculus_40MD33", parallel.core = parallel::detectCores() - 1,
                                           verbose = TRUE)
##
##
############################################
############     Mugil liza     ############
############################################
##
##Keeping the population
liza33 <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("liza_SA"))
##
##Filtering
##
#Filter on reproducibility, threshold (here called t, do not ask why) 97% reproducible - a meassurement of quality per loci
liza33_rep<-gl.filter.repavg(liza33, t=0.97)
liza33_rep
#Fitler by monomorphic loci - as they do not provide information for population structure and simply slow the analysis
liza33_mono<-gl.filter.monomorphs(liza33_rep, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
liza33_mono
#Filter by paralogs
liza33_trim<- gl.filter.hamming(liza33_mono, threshold=0.2, pb=T, v=5)
liza33_trim
##################
#Filter 0% missing data
liza33_0MD <- gl.filter.callrate(liza33_trim, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
liza33_0MD
#Filter to remove all but one of multiple snps in the same fragment
liza33_0MD_1SNP<-gl.filter.secondaries(liza33_0MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
liza33_0MD_1SNP
##
##Convert to VCF using radiator package
##
library(radiator)
vcf_liza33_0MD<-genomic_converter(liza33_0MD_1SNP, strata = NULL, output = "vcf",
                                        filename = "liza_0MD33", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
##################
#Filter 40% missing data
liza33_40MD <-gl.filter.callrate(liza33_trim, method = "loc", threshold = 0.6, mono.rm = F, recalc = F, plot = F, v = 2)
liza33_40MD
#Filter to remove all but one of multiple snps in the same fragment
liza33_40MD_1SNP<-gl.filter.secondaries(liza33_40MD, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
liza33_40MD_1SNP
##Convert to VCF using radiator package
##
vcf_liza33_40MD<-genomic_converter(liza33_40MD_1SNP, strata = NULL, output = "vcf",
                                           filename = "liza_40MD33", parallel.core = parallel::detectCores() - 1,
                                           verbose = TRUE)
##
##
##Convert the VCF files into DaDi’s format (.snp) using the python script available at https://github.com/CoBiG2/RAD_Tools/blob/master/vcf2DaDi.py. 
#
#################################################################################################################################################
#
#jessika.neves@icbs.ufal.br