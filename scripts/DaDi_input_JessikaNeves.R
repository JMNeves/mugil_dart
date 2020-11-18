########################
###                  ###
###   DartR - DADI   ###
### By Jessika Neves ###
###   November 2020  ###
###                  ###
########################
#
#This script was developed under the dartR version 1.1.6 available at https://cran.r-project.org/src/contrib/Archive/dartR/
#
##Producing a dataset to run DaDi
##As the analysis requeires one dataset per species, it is needed to filter by species before the other filters
#
library(dartR)
library(radiator)
##
############################################
############ Mugil brevirostris ############
############################################
##
##Keeping the population
brevirostris <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("brev_SA"))
##
##Filtering
##
#Filter 0% missing data
brevirostris_0MD <- gl.filter.callrate(brevirostris, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility
brevirostris_rep<-gl.filter.repavg(brevirostris_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
brevirostris_snp<-gl.filter.secondaries(brevirostris_rep, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
#Fitler by monomorphic loci 
brevirostris_mono<-gl.filter.monomorphs(brevirostris_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
brevirostris_hwe<-gl.filter.hwe(brevirostris_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
brevirostris_trim<- gl.filter.hamming(brevirostris_hwe, threshold=0.2, pb=T, v=5)
##
##Convert to VCF using radiator package
##
vcf_brevirostris_0MD<-genomic_converter(brevirostris_trim, strata = NULL, output = "vcf",
                                        filename = "brevirostris_0MD", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
vcf_brevirostris_0MD$tidy.data
tidy_brevirostris_0MD<-vcf_brevirostris_0MD$tidy.data
write.matrix(tidy_brevirostris_0MD,"~/dartr/brevirostris/tidy_brevirostris_0MD.txt")
##
##
############################################
############    Mugil curema    ############
############################################
##
##Keeping the population
curema <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("cure_SA", "inci_SA"))
##
##Filtering
##
#Filter 0% missing data
curema_0MD <- gl.filter.callrate(curema, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility
curema_rep<-gl.filter.repavg(curema_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
curema_snp<-gl.filter.secondaries(curema_rep, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
#Fitler by monomorphic loci 
curema_mono<-gl.filter.monomorphs(curema_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
curema_hwe<-gl.filter.hwe(curema_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
curema_trim<- gl.filter.hamming(curema_hwe, threshold=0.2, pb=T, v=5)
##
##Convert to VCF using radiator package
##
vcf_curema_0MD<-genomic_converter(curema_trim, strata = NULL, output = "vcf",
                                        filename = "curema_0MD", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
vcf_curema_0MD$tidy.data
tidy_curema_0MD<-vcf_curema_0MD$tidy.data
write.matrix(tidy_curema_0MD,"~/dartr/curema/tidy_curema_0MD.txt")
##
##
############################################
############   Mugil curvidens  ############
############################################
##
##Keeping the population
curvidens <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("curv_SA", "curv_MB"))
##
##Filtering
##
#Filter 0% missing data
curvidens_0MD <- gl.filter.callrate(curvidens, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility
curvidens_rep<-gl.filter.repavg(curvidens_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
curvidens_snp<-gl.filter.secondaries(curvidens_rep, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
#Fitler by monomorphic loci 
curvidens_mono<-gl.filter.monomorphs(curvidens_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
curvidens_hwe<-gl.filter.hwe(curvidens_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
curvidens_trim<- gl.filter.hamming(curvidens_hwe, threshold=0.2, pb=T, v=5)
##
##Convert to VCF using radiator package
##
vcf_curvidens_0MD<-genomic_converter(curvidens_trim, strata = NULL, output = "vcf",
                                        filename = "curvidens_0MD", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
vcf_curvidens_0MD$tidy.data
tidy_curvidens_0MD<-vcf_curvidens_0MD$tidy.data
write.matrix(tidy_curvidens_0MD,"~/dartr/curvidens/tidy_curvidens_0MD.txt")
##
##
############################################
############  Mugil rubrioculus ############
############################################
##
##Keeping the population
rubrioculus <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("rubr_SA"))
##
##Filtering
##
#Filter 0% missing data
rubrioculus_0MD <- gl.filter.callrate(rubrioculus, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility
rubrioculus_rep<-gl.filter.repavg(rubrioculus_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
rubrioculus_snp<-gl.filter.secondaries(rubrioculus_rep, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
#Fitler by monomorphic loci 
rubrioculus_mono<-gl.filter.monomorphs(rubrioculus_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
rubrioculus_hwe<-gl.filter.hwe(rubrioculus_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
rubrioculus_trim<- gl.filter.hamming(rubrioculus_hwe, threshold=0.2, pb=T, v=5)
##
##Convert to VCF using radiator package
##
vcf_rubrioculus_0MD<-genomic_converter(rubrioculus_trim, strata = NULL, output = "vcf",
                                        filename = "rubrioculus_0MD", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
vcf_rubrioculus_0MD$tidy.data
tidy_rubrioculus_0MD<-vcf_rubrioculus_0MD$tidy.data
write.matrix(tidy_rubrioculus_0MD,"~/dartr/rubrioculus/tidy_rubrioculus_0MD.txt")
##
##
############################################
############     Mugil liza     ############
############################################
##
##Keeping the population
liza <- gl.keep.pop(dart, mono.rm = FALSE, pop.list=c("liza_SA"))
##
##Filtering
##
#Filter 0% missing data
liza_0MD <- gl.filter.callrate(liza, method = "loc", threshold = 1.0, mono.rm = F, recalc = F, plot = F, v = 2)
#Filter on reproducibility
liza_rep<-gl.filter.repavg(liza_0MD, t=0.97)
#Filter to remove all but one of multiple snps in the same fragment
liza_snp<-gl.filter.secondaries(liza_rep, method="best") #filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order (method="best") or at random (method="random")
#Fitler by monomorphic loci 
liza_mono<-gl.filter.monomorphs(liza_snp, v=5) #v=verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#Filter by HWE
liza_hwe<-gl.filter.hwe(liza_mono, alpha = 0.05, basis = "any", bon = TRUE, v = 5)
#Filter by paralogs
liza_trim<- gl.filter.hamming(liza_hwe, threshold=0.2, pb=T, v=5)
##
##Convert to VCF using radiator package
##
vcf_liza_0MD<-genomic_converter(liza_trim, strata = NULL, output = "vcf",
                                        filename = "liza_0MD", parallel.core = parallel::detectCores() - 1,
                                        verbose = TRUE)
vcf_liza_0MD$tidy.data
tidy_liza_0MD<-vcf_liza_0MD$tidy.data
write.matrix(tidy_liza_0MD,"~/dartr/liza/tidy_liza_0MD.txt")
##
##
##Convert the files into DaDi’s format (.snp) using the python script available at https://github.com/CoBiG2/RAD_Tools/blob/master/vcf2DaDi.py. 
#
#################################################################################################################################################
#
#jessika.neves@icbs.ufal.br
