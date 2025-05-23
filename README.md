# CompMap
Compares reads mapped to two references and counts reads for allele-specific expression

## Dependencies

1. Python 3.x
2. [pysam](https://pysam.readthedocs.io/en/latest/api.html)
3. art (for pretty ASCII art)

## Installation

### Direct and easy installation with `pip`

`pip` will install all the dependencies and make a path-ready excecutable.

    python3 -m pip install https://github.com/santiagosnchez/CompMap/raw/master/dist/CompMap-1.1-py3-none-any.whl

also:

    python3 -m pip install --user https://github.com/santiagosnchez/CompMap/raw/master/dist/CompMap-1.1.tar.gz

works if your systems does not take `wheel` files. Use the `--user` flag if you want to install locally.

### Manual installation and download

The program requires `pysam` to be installed. To do that just run `pip`:

    pip install pysam

To install locally:

    pip install --user pysam

Other dependencies can be installed in the same way:

    pip install art time

Tu get just the code, run:

    wget https://raw.githubusercontent.com/santiagosnchez/CompMap/master/CompMap

For the whole repo:

    git clone https://github.com/santiagosnchez/CompMap.git

## Workflow and rationale

**Note:** Currently, `CompMap` does not accept the union of multiple features such as exons, **nor will it work with paired-ended reads** (future work).

The basic workflow of the program consists of parsing two BED files (one for each reference to which the mixed-read data was aligned to) with the coordinates of the features of interest, such as genes or transcripts. Two BAM files will be parsed: each one aligned to two different references (i.e. two different species or genotypes). Reads with better mapping features (namely 'alignment score' and 'number of mismatches') in a given BAM file will be counted towards a given reference. However, only a fraction `prop` of ambiguously aligning reads (i.e. with the same 'alignment score' and/or 'number of mismatches') will be counted.

Generally, `ambiguous_reads_in_sp1 == ambiguous_reads_in_sp2`, so for now we can call them:

```
ambiguous_matches
```

But once we have counts for the number of best matches and the total count that includes ambiguous reads in gene x as:

```
best_matches_sp1
best_matches_sp2
total_matches_sp1
total_matches_sp2
```

You can calculate the fraction of `ambiguous_matches` that contributes to genotype "sp1" as:

```
if best_matches_sp1 < best_matches_sp2:
    # get the proportional difference between both counts
    prop = best_matches_sp1/best_matches_sp2
    # add the surplus counts in the proportional difference
    corrected_sp1 = best_matches_sp1 + ((total_matches_sp1-best_matches_sp1) * prop)
    corrected_sp2 =  best_matches_sp2 + ((total_matches_sp2-best_matches_sp2) * 1-prop)
if best_matches_sp1 > best_matches_sp2
    # do the same both with:
    prop = best_matches_sp2/best_matches_sp1
    corrected_sp1 = best_matches_sp1 + ((total_matches_sp1-best_matches_sp1) * 1-prop)
    corrected_sp2 =  best_matches_sp2 + ((total_matches_sp2-best_matches_sp2) * prop)
```

The alignment stats that `CompMap` looks at are:

1. Alignment score ('AS')
2. Number of mismatches ('nM')

But BAM-specific \"alignment score\" and \"number of mismatches\" tags can be specified as arguments.

## Download/installation

For just the code:

    wget https://raw.githubusercontent.com/santiagosnchez/CompMap/master/CompMap.py

For the whole repo:

    git clone https://github.com/santiagosnchez/CompMap.git

## Arguments

Run the program with the `-h` or `-help` argument to look at the different options and examples.

```
  ____                           __  __               
 / ___|  ___   _ __ ___   _ __  |  \/  |  __ _  _ __  
| |     / _ \ | '_ ` _ \ | '_ \ | |\/| | / _` || '_ \
| |___ | (_) || | | | | || |_) || |  | || (_| || |_) |
 \____| \___/ |_| |_| |_|| .__/ |_|  |_| \__,_|| .__/
                         |_|                   |_|    

usage: CompMap [-h] --bam1 BAM1 --bam2 BAM2 --bed1 BED1 --bed2 BED2
               [--base BASE] [--AS_tag AS_TAG] [--NM_tag NM_TAG] [--star]
               [--binom] [--suffix1 SUFFIX1] [--suffix2 SUFFIX2]

    Compares read mapping features of two BAM files  aligned to two different references
    and generates reference-specific read counts based on best matches.

    This program is useful for counting reads from mixed-read data such as allele-specific expression RNA-seq assays.

optional arguments:
  -h, --help            show this help message and exit
  --bam1 BAM1, -1 BAM1  1st bam file. Preferentially sorted and indexed.
  --bam2 BAM2, -2 BAM2  2nd bam file. Preferentially sorted and indexed.
  --bed1 BED1, -b1 BED1
                        first bed file. The 4th column is expected to be the gene id alone and must be the same as in bed2.
  --bed2 BED2, -b2 BED2
                        second bed file. The 4th column is expected to be the gene id alone and must be the same as in bed1.
  --base BASE, -b BASE  base name for your output files (recommended).
  --AS_tag AS_TAG       provide an "alignment score" tag in your BAM file. Decault: AS
  --NM_tag NM_TAG       provide a "number of mismatches tag" tag in your BAM file. Default: nM
  --star                use the NH tag to filter out multi-mapping reads.
  --binom               use binomial distribution to sample ambiguous ASE read counts.
  --suffix1 SUFFIX1, -s1 SUFFIX1
                        string suffix used to distinguish matches, mismatches, and ambiguous reads in -1. Default: 1
  --suffix2 SUFFIX2, -s2 SUFFIX2
                        string suffix used to distinguish matches, mismatches, and ambiguous reads in -2. Default: 2

Examples:
CompMap -1 species1.bam -2 species2.bam -b1 genes_sp1.bed -b2 genes_sp2.bed -b sample_1 -s1 sp1 -s2 sp2 --NM_tag NM

Example of BED file:
chr1	0	1000	gene1
chr1	2500	3000	gene2
chr2	0	2500	gene3
chr3	150	1000	gene4
```

## Test data

The `test_data` directory in this repo includes the files necessary to run a test with 100 simulated genes and RNA-seq data.

To run the test use the following code:

```
cd test_data
../CompMap -1 sample_01_d1.short.bam \
           -2 sample_01_d2.short.bam \
           -b1 highrate_d0.1_100genes_d1.bed \
           -b2 highrate_d0.1_100genes_d2.bed \
           --base sample_01 \
           --NM_tag NM \
           -s1 d1 \
           -s2 d2
```

You should end up with a text file named `sample_01_counts.txt` that includes allele-specifc counts.

## Running simulations

The repo also includes code to simulate your own data set. The scripts can be found under the `cds_and_rnaseq_simulations` directory. The simulation pipeline can be found in the `1_pipeline_sim_reads.sh` script. The requirements to run this simulation pipeline properly are:

* gnu-parallel
* samtools
* bwa
* featureCounts
* CompMap

Once count files are generated, the `R` script `DE_ASE_simulated_data.R` includes code to analyze the data using `DESeq2`, `tidyr`, and `ggplot2` for visualization.
