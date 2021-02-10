#!/usr/local/bin/Rscript
library(polyester)
library(Biostrings)
library(ggplot2)
library(tidyr)

# read command arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) !=2 ){
    stop("Fasta files for both sets of transcripts must be specified.", call.=FALSE)
}

# make a function to return a matrix with col = number of treatments
# based on a vector with directional fold-changes
# if fc[1] == 5, then fold_changes[1,] == c(6, 1)
# if fc[1] == -5, then fold_changes[1,] == c(1, 6)
# if fc[1] == 0, then fold_changes[1,] == c(1, 1)
fc2matrix = function(fc, ncol=2, nrow=1000){
    fold_changes = matrix(1, ncol=ncol, nrow=nrow)
    fold_changes[fc > 0, 2] = 1
    fold_changes[fc > 0, 1] = fc[ fc > 0 ] + 1
    fold_changes[fc < 0, 1] = 1
    fold_changes[fc < 0, 2] = abs(fc[ fc < 0 ]) + 1
    return(fold_changes)
}

#########################################
# expression divergence between species #
#########################################

# keep seed for consistency between simulations
set.seed(1)
# sample from an exponential distribution with a rate = 0.5
# most values will be close to 0 and a few will be divergent in expression
species_fc = round(rexp(1000, 0.5),0)
# decide the direction in fold-change based on a binomial distribution with prob = 0.5
# keep seed for consistency between simulations
set.seed(1)
species_fc_neg = rbinom(1000, 1, 0.5) * -1
species_fc_neg[ species_fc_neg == 0 ] = 1
# vector with directional fold-changes
species_fc = species_fc * species_fc_neg
# matrix with fold-changes
species_fold_changes = fc2matrix(species_fc)

##############################
# allele-specific expression #
##############################

# keep seed for consistency between simulations
set.seed(2)
# sample from an exponential distribution with a rate = 0.5
# most values will be close to 0 and a few will be divergent in expression
allele_fc = round(rexp(1000, 0.7),0)
# decide the direction in fold-change based on a binomial distribution with prob = 0.5
# keep seed for consistency between simulations
set.seed(2)
allele_fc_neg = rbinom(1000, 1, 0.5) * -1
allele_fc_neg[ allele_fc_neg == 0 ] = 1
# vector with directional fold-changes
allele_fc = allele_fc * allele_fc_neg
# check distribution of fold-changes
hist(allele_fc)
# matrix with fold-changes
allele_fold_changes = fc2matrix(allele_fc)

####################################
# combine fold-change and classify #
####################################

fold_change_comb = as.data.frame(cbind(species_fc, allele_fc))
real_class = rep(NA, 1000)
real_class[ fold_change_comb[,1] == 0 & fold_change_comb[,2] == 0 ] = "conserved"
real_class[ is.na(real_class) & fold_change_comb[,1] == fold_change_comb[,2] ] = "cis only"
real_class[ is.na(real_class) & (fold_change_comb[,1] != 0 & fold_change_comb[,2] == 0) ] = "trans only"
real_class[ is.na(real_class) & (fold_change_comb[,1] == 0 & fold_change_comb[,2] != 0) ] = "cis-trans (compensatory)"
real_class[ is.na(real_class) & ((fold_change_comb[,1] > 0 & fold_change_comb[,2] < 0) | (fold_change_comb[,1] < 0 & fold_change_comb[,2] > 0)) ] = "cis-trans x (compensatory)"
real_class[ is.na(real_class) & ((fold_change_comb[,1] > 0 & fold_change_comb[,2] > 0) | (fold_change_comb[,1] < 0 & fold_change_comb[,2] < 0)) ] = "cis-trans (enhancing)"
fold_change_comb = cbind(fold_change_comb, real_class)

# save data:
#rownames(fold_change_comb) = names(briggsae_fasta_small)
write.csv(fold_change_comb, file="sim_fold_change_comb.csv", quote=FALSE)

# melt df
fold_change_comb_melt = gather(fold_change_comb, species_vs_allele, fc, -real_class)
# plot and check fold-changes
dev.new()
ggplot(fold_change_comb_melt, aes(x=fc)) + geom_histogram(fill="darkgreen", alpha=0.6) + facet_wrap(~species_vs_allele)

# plot real classes
dev.new()
ggplot(fold_change_comb, aes(x=species_fc, y=allele_fc, color=real_class)) + geom_point(alpha=0.5) + geom_jitter()

##########################
## Process sequence data #
##########################

# CDS sequence fasta files
d1_fasta_file = args[1]
d2_fasta_file = args[2]

# read objects
d1_fasta = readDNAStringSet(d1_fasta_file)
d2_fasta = readDNAStringSet(d2_fasta_file)

# sample 1000 genes
# set.seed(1)
# idx = sample(1:length(briggsae_fasta), 1000)
# briggsae_fasta_small = briggsae_fasta[idx]
# nigoni_fasta_small = nigoni_fasta[idx]

# save to file
# writeXStringSet(briggsae_fasta_small, 'briggsae.cds.small.fa')
# writeXStringSet(nigoni_fasta_small, 'nigoni.cds.small.fa')

##################################
## simulate between species data #
##################################

# baseline number of reads with 40x coverage
readspertx = round(20 * width(d1_fasta) / 100)

# output
outdir_1 = sub("\\.fasta","", d1_fasta_file)
outdir_2 = sub("\\.fasta","", d2_fasta_file)

# simulate briggsae reads
# three replicates each
print("simulating between species data")
simulate_experiment(d1_fasta_file,
                   reads_per_transcript=readspertx,
                   #meanmodel=TRUE,
                   num_reps=c(3,3),
                   fold_changes=species_fold_changes,
                   outdir=paste0('species_',outdir_1,'_simulated_reads'),
                   paired=FALSE)
print(paste("done for",d1_fasta_file))
# simulate nigoni reads
simulate_experiment(d2_fasta_file,
                   reads_per_transcript=readspertx,
                   num_reps=c(3,3),
                  #meanmodel=TRUE,
                   fold_changes=species_fold_changes,
                   outdir=paste0('species_',outdir_2,'_simulated_reads'),
                   paired=FALSE)
print(paste("done for",d2_fasta_file))

##################################
## simulate between alleles data #
##################################

# baseline number of reads with 20x coverage
# allele_readspertx_briggsae = round(20 * width(briggsae_fasta_small) / 100)
# allele_readspertx_nigoni = round(20 * width(nigoni_fasta_small) / 100)

# simulate briggsae reads
# three replicates each
print("simulating between allele data")
simulate_experiment(d1_fasta_file,
                   reads_per_transcript=readspertx,
                   num_reps=c(3,3),
                   #meanmodel=TRUE,
                   fold_changes=allele_fold_changes,
                   outdir=paste0('alleles_',outdir_1,'_simulated_reads'),
                   paired=FALSE)
print(paste("done for",d1_fasta_file))
# simulate nigoni reads
simulate_experiment(d2_fasta_file,
                   reads_per_transcript=readspertx,
                   num_reps=c(3,3),
                   #meanmodel=TRUE,
                   fold_changes=allele_fold_changes,
                   outdir=paste0('alleles_',outdir_2,'_simulated_reads'),
                   paired=FALSE)
print(paste("done for",d2_fasta_file))
