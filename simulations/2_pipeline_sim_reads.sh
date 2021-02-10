#!/bin/bash

#################################
# simulate diverged CDS regions #
#################################

Rscript run_cds_sim.R

###################
# run simulations #
###################

# simulate RNA-seq data
ls *.fasta | xargs -n 2 Rscript simulate_rnaseq.R

# compress fasta files
#find . -type f -name "*.fasta" | grep simulated | parallel gzip {}

##########################
# Process simulated data #
##########################

# make unique directories for each simulation
ls -d *simulated_reads | sed -E 's/d[12]_|comb_//' | uniq | xargs mkdir

# delete same allele reads
ls *d1_simulated_reads/*.fasta | grep -E "sample_0[1-3]" | parallel --colsep="/" --plus mv {1}/{2} {1/d1_/}/{2.}_d1.fasta
ls *d2_simulated_reads/*.fasta | grep -E "sample_0[4-6]" | parallel --colsep="/" --plus mv {1}/{2} {1/d2_/}/{2.}_d2.fasta
# ls *d1_simulated_reads/*.bam | grep -E "sample_0[1-3]" | parallel --colsep="/" --plus mv {1}/{2} {1/d1_/}/{2.}_d1.bam
# ls *d2_simulated_reads/*.bam | grep -E "sample_0[4-6]" | parallel --colsep="/" --plus mv {1}/{2} {1/d2_/}/{2.}_d2.bam

# delete old dirs
ls -d *_simulated_reads | grep -E "d1|d2" | xargs rm -rf

# map to reference
find . -type f -name "*.fasta" | grep simulated | \
cut -d"/" -f2,3 | tr '/' ' ' | sed -E 's/_simulated_reads//' | \
parallel --colsep=" " --plus "bwa mem bwa_index/{1/species_|alleles_/}.fasta {1}_simulated_reads/{2} | samtools view -b - > {1}_simulated_reads/{2.}.bam"

# combine alleles
ls -d alleles* | sed 's/d[12]/comb/' | uniq | xargs mkdir
ls -d alleles* | xargs -n 3 | parallel --colsep=" " "cat {2}/sample_01.fasta {3}/sample_04.fasta > {1}/sample_01.fasta"
ls -d alleles* | xargs -n 3 | parallel --colsep=" " "cat {2}/sample_02.fasta {3}/sample_05.fasta > {1}/sample_02.fasta"
ls -d alleles* | xargs -n 3 | parallel --colsep=" " "cat {2}/sample_03.fasta {3}/sample_06.fasta > {1}/sample_03.fasta"
# make subdirectories for each combined simulation and each dataset
# ls -d alleles_*comb* | xargs -I% sh -c "mkdir %/d1 && mkdir %/d2"

# combined data: map to reference
ls alleles_*comb*/*.fasta | cut -d"/" -f1,2 | tr '/' ' ' | sed -E 's/_comb_simulated_reads//' | \
parallel --colsep=" " --plus "bwa mem bwa_index/{1/alleles_/}_d1.fasta {1}_comb_simulated_reads/{2} | samtools view -b - > {1}_comb_simulated_reads/{2.}_d1.bam"
ls alleles_*comb*/*.fasta | cut -d"/" -f1,2 | tr '/' ' ' | sed -E 's/_comb_simulated_reads//' | \
parallel --colsep=" " --plus "bwa mem bwa_index/{1/alleles_/}_d2.fasta {1}_comb_simulated_reads/{2} | samtools view -b - > {1}_comb_simulated_reads/{2.}_d2.bam"

# remove fasta files
find . -type f -name "*.fasta" | grep simulated | xargs rm

# counts with featureCounts
mkdir counts
find . -type f -name "*.bam" | grep -v comb | sort | xargs -n 3 | \
sed -E 's/\.\/|_simulated_reads.*_(d[12])\.bam.*/_\1/g' | sed 's/_//' | xargs -I% echo % % | sed -E 's/alleles_|species_//' > gffs
find . -type f -name "*.bam" | grep -v comb | sort | xargs -n 3 | paste -d' ' gffs - | \
parallel --colsep=" " featureCounts -a {1}.gff -F GFF -t transcript -o counts/{2}_counts.txt {3} {4} {5} &> /dev/null &
rm gffs

# counts with CompMap
# sort and index BAM files
ls alleles_*comb*/*.bam | parallel samtools sort -o {.}.sort.bam {}
ls alleles_*comb*/*.bam | xargs -n 2 | parallel --colsep=" " mv {2} {1}
ls alleles_*comb*/*.bam | parallel samtools index {}

# prepare argument list
ls alleles_*comb*/*.bam | xargs -n 2 > tmp1
cat tmp1 | sed -E 's/alleles_|_comb_.*_0[1-3]|\.bam//g' | sed 's/$/\.bed/' | paste -d' ' - | xargs -I% echo % % | sed 's/d2/d1/' | paste -d' ' - tmp1 > tmp2
cut -d' ' -f 1,3 tmp2 | sed 's/d[12]\.bed//' | sed -E 's/ alleles.*\/|_d1\.bam//g' | paste -d' ' - tmp2 > tmp3
cat tmp1 | grep -o -E "d[12]" | xargs -n 2 | paste -d' ' - tmp3 > tmp4

#ls *.bed | xargs -I % sh -c "cut -f1 highrate_1000_0.1_d1.bed | sed 's/d1\.//' | paste % - > tmp && mv tmp %"

# run CompMap
cat tmp4 | parallel --colsep=" " CompMap -s1 {1} -s2 {2} --base {3} -b1 {4} -b2 {5} --NM_tag NM  -1 {6} -2 {7}
rm tmp*
mv *_counts.txt counts

# prepare/consolidate count files
cd counts
mkdir final
paste species_highrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|highrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/species_highrate_counts.txt
paste species_mediumrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|mediumrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/species_mediumrate_counts.txt
paste species_lowrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|lowrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/species_lowrate_counts.txt
paste alleles_highrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|highrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/alleles_highrate_counts.txt
paste alleles_mediumrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|mediumrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/alleles_mediumrate_counts.txt
paste alleles_lowrate_*_d[12]*.txt | tail -n +2 | perl -pe 's/\.\/|lowrate.+?\/|\.bam//g' | cut -f1,7-9,16-18 > final/alleles_lowrate_counts.txt

# consolidate CompMap counts
paste highrate_*_sample_0[1-3]_* | perl -pe 's/highrate.+?_sample/sample/g' | cut -f1-3,5-6,8-9 > final/compmap_highrate_counts.txt
paste mediumrate_*_sample_0[1-3]_* | perl -pe 's/mediumrate.+?_sample/sample/g' | cut -f1-3,5-6,8-9 > final/compmap_mediumrate_counts.txt
paste lowrate_*_sample_0[1-3]_* | perl -pe 's/lowrate.+?_sample/sample/g' | cut -f1-3,5-6,8-9 > final/compmap_lowrate_counts.txt

# archive old files
mkdir archive
mv *txt* archive
cd ../
# create a directory for figures
mkdir figures

# end of pipeline
