#!python env
import sys
import argparse
#import itertools
import art
import time
try:
    import pysam
except ImportError as er:
    print("""please install pysam:
pip install pysam\n""")
    print(er)
    sys.exit()

def main():
    # parse arguments
    pretty_title = art.text2art("CompMap")
    print(pretty_title)
    parser = argparse.ArgumentParser(prog="CompMap.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    Compares reads from one sample aligned to two different references
    and sorts out the best match for each read to two distinct bam files.

    This program is useful for splitting reads from mixed-read data (e.g. from hybrids).\n""",
    epilog="""
Examples:
python CompMap.py -1 first.bam -2 second.bam -r read_name_list -b my_base_name -s1 sp1 -s2 sp2

Example of BED file:
chr1\t0\t1000\tgene1
chr1\t2500\t3000\tgene2
chr2\t0\t2500\tgene3
chr3\t150\t1000\tgene4
    \n""")
    parser.add_argument(
    '--bam1', '-1', type=str, required=True,
    help='1st bam file. Preferentially sorted and indexed.')
    parser.add_argument(
    '--bam2', '-2', type=str, required=True,
    help='2nd bam file. Preferentially sorted and indexed.')
    parser.add_argument(
    '--bed1', '-b1', type=str, required=True,
    help='first bed file. The 4th column is expected to be the gene id alone and must be the same as in bed2.')
    parser.add_argument(
    '--bed2', '-b2', type=str, required=True,
    help='second bed file. The 4th column is expected to be the gene id alone and must be the same as in bed1.')
    parser.add_argument(
    '--base', '-b', type=str, default="out",
    help='base name for your output files (recommended).')
    parser.add_argument(
    '--AS_tag', type=str, default="AS",
    help='provide an \"alignment score\" tag in your BAM file. Decault: AS')
    parser.add_argument(
    '--NM_tag', type=str, default="nM",
    help='provide a \"number of mismatches tag\" tag in your BAM file. Default: nM')
    parser.add_argument(
    '--suffix1', '-s1', type=str, default="_1",
    help='string suffix used to distinguish matches, mismatches, and ambiguous reads in -1. Default: 1')
    parser.add_argument(
    '--suffix2', '-s2', type=str, default="_2",
    help='string suffix used to distinguish matches, mismatches, and ambiguous reads in -2. Default: 2')
    args = parser.parse_args()

    # read bam files
    print("Reading BAM files.")
    bam1 = pysam.AlignmentFile(args.bam1, 'rb')
    bam2 = pysam.AlignmentFile(args.bam2, 'rb')

    # index by read name
    print("Indexing file by read name:",args.bam1)
    name_indexed1 = index_by_readName(bam1)
    print("Indexing file by read name:",args.bam2)
    name_indexed2 = index_by_readName(bam2)

    # copy header
    header1 = bam1.header.copy()
    header2 = bam2.header.copy()
    # output files
    out = (args.base,args.suffix1,args.suffix2)

    # read bed files
    print("Reading BED files")
    bed1 = read_bed(args.bed1)
    bed2 = read_bed(args.bed2)

    # compare read counts, main function
    comp_map(bed1, bed2, bam1, bam2, name_indexed1, name_indexed2, out, args.AS_tag, args.NM_tag)

# function to compare reads
def comp_map(bed1, bed2, bam1, bam2, name_indexed1, name_indexed2, out, as_tag, nm_tag):
    # process bed files
    # check if gene ids are the same
    print("Checking gene_id names: ", end="")
    bed1_gn = [ x[3] for x in bed1 ]
    bed2_gn = [ x[3] for x in bed2 ]
    bed1_gn.sort()
    bed2_gn.sort()
    if bed1_gn != bed2_gn:
        raise NameError('Gene ids in BED files are not the same.')
        sys.exit(1)
    print("pass.")
    total_features = len(bed1_gn)

    # change to dictionaries
    bed1 = { x[3]:x for x in bed1 }
    bed2 = { x[3]:x for x in bed2 }

    # start counter
    feature_counter = 0

    # open outfile and write header
    print("Writting to: "+out[0]+"_counts.txt")
    o = open(out[0]+"_counts.txt", "w")
    o.write("gene_id\t"+out[0]+"_"+out[1]+"\t"+out[0]+"_"+out[2]+"\n")

    # loop through bedfiles
    t1 = time.time()
    for geneid in bed1_gn:
        feature_counter += 1
        bedline1 = bed1[geneid]
        bedline2 = bed2[geneid]
        # count reads
        counts1 = cm_count_reads(bedline1, bam1, name_indexed1, bam2, name_indexed2, as_tag, nm_tag)
        counts2 = cm_count_reads(bedline2, bam2, name_indexed2, bam1, name_indexed1, as_tag, nm_tag)
        # correct counts
        corr1, corr2 = correct_counts(counts1, counts2)
        o.write(geneid+"\t"+str(corr1)+"\t"+str(corr2)+"\n")

        # print progress bar
        progress = make_progress_bar(feature_counter, total_features, t1, 70)
        print("\r", progress[0] % progress[1:], end='', flush=True)
        # if feature_counter % 100 == 0:
        #     print("processed "+str(feature_counter)+" features")
    o.close()
    print("\nDone processing "+str(feature_counter)+" features.")

def cm_count_reads(bedline, bam, bam_ni, bam_alt, bam_ni_alt, as_tag, nm_tag):
    # check unique alignment
    def check_unique(read_name, bam_ni):
        c = 0
        for r in bam_ni.find(read_name):
            c += 1
        if c > 1:
            return False
        else:
            return True
    # check if read is present
    def check_is_found(read_name, bam_ni):
        try:
            bam_ni.find(read_name)
        except KeyError:
            return False
        else:
            return True
    # return score and mismatches
    def get_as_nm(x, as_tag, nm_tag):
        return x.get_tag(as_tag), x.get_tag(nm_tag)
    # compare reads based on AS and NM
    def compare_reads(x, x_alt, as_tag, nm_tag):
        score, mismatches = get_as_nm(x, as_tag, nm_tag)
        score_alt, mismatches_alt = get_as_nm(x_alt, as_tag, nm_tag)
        if score > score_alt:
            return 'M'
        elif mismatches < mismatches_alt:
            return 'M'
        else:
            if score == score_alt:
                return 'A'
            elif mismatches == mismatches_alt:
                return 'A'
    #
    count_bmt = 0; count_amb = 0
    for r in bam.fetch(bedline[0], int(bedline[1]), int(bedline[2])):
        if not any([r.is_unmapped, r.is_secondary, r.is_supplementary]) and check_unique(r.qname, bam_ni):
            if check_is_found(r.qname, bam_ni_alt) and check_unique(r.qname, bam_ni_alt):
                r_alt = next(bam_ni_alt.find(r.qname))
                if not any([r_alt.is_unmapped, r_alt.is_secondary, r_alt.is_supplementary]):
                    cr = compare_reads(r, r_alt, as_tag, nm_tag)
                    if cr == 'M':
                        count_bmt += 1
                    elif cr == 'A':
                        count_amb += 1
                else:
                    count_bmt += 1
            else:
                count_bmt += 1
    count_all = count_bmt + count_amb
    return count_all, count_bmt, count_amb

def correct_counts(counts1, counts2):
    if counts1[0] == counts1[1] and counts2[0] == counts2[1]:
        return counts1[0], counts2[0]
    elif counts1[0] == counts1[1] and counts2[0] > counts2[1]:
        return counts1[0], counts2[0]+counts2[1]
    elif counts1[0] > counts1[1] and counts2[0] == counts2[1]:
        return counts1[0]+counts1[1], counts2[0]
    else:
        if counts1[1] < counts2[1]:
            prop = counts1[1]/counts2[1]
            corrected1 = counts1[1] + ((counts1[0]-counts1[1]) * prop)
            corrected2 = counts2[1] + ((counts2[0]-counts2[1]) * 1-prop)
            return int(corrected1), int(corrected2)
        elif counts1[1] > counts2[1]:
            prop = counts2[1]/counts1[1]
            corrected1 = counts1[1] + ((counts1[0]-counts1[1]) * 1-prop)
            corrected2 = counts2[1] + ((counts2[0]-counts2[1]) * prop)
            return int(corrected1), int(corrected2)
        else:
            corrected1 = counts1[1] + ((counts1[0]-counts1[1]) * .5)
            corrected2 = counts2[1] + ((counts2[0]-counts2[1]) * .5)
            return int(corrected1), int(corrected2)

# index bam by read name
def index_by_readName(bam):
    name_indexed = pysam.IndexedReads(bam)
    name_indexed.build()
    return name_indexed

def read_bed(file):
    with open(file, "r") as f:
        bed = [ x.split("\t") for x in f.read().splitlines() ]
    return bed

def make_progress_bar(rec, total, t1, width):
    i = (rec/total * 100) % 100
    if i != 0:
        asterix = "*" * int(i * (width/100))
        dashes = "-" * (width-int(i * width/100))
    else:
        asterix = "*" * width
        dashes = ""
        i = 100
    t2 = time.time()
    elapsed = t2-t1
    return "|"+asterix+dashes+"| "+"%5.2f%% %7.2f s", i, elapsed

if __name__ == "__main__":
    main()