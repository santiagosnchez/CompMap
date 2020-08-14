#!/usr/bin/python env
import sys
import argparse
try:
    import pysam
except ImportError as er:
    print("""please install pysam:
pip install pysam\n""")
    print(er)
    sys.exit()

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="CompMap.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    Compares reads from one sample aligned to two different references
    and sorts out the best match for each read to two distinct bam files.

    This program is useful for splitting reads from mixed-read data (e.g. from hybrids).\n""",
    epilog="""
Examples:
python CompMap.py -1 first.bam -2 second.bam -r read_name_list -b my_base_name -s1 sp1 -s2 sp2

Example for read list:
ABC-HG000:000:XXXXXXX:1:0010:001:100
ABC-HG000:000:XXXXXXX:1:0010:001:110
ABC-HG000:000:XXXXXXX:1:0010:001:130
ABC-HG000:000:XXXXXXX:1:0010:001:110
    \n""")
    parser.add_argument(
    '--bam1', '-1', type=str, required=True,
    help='first bam file.')
    parser.add_argument(
    '--bam2', '-2', type=str, required=True,
    help='second bam file.')
    parser.add_argument(
    '--reads', '-r', type=str, required=True,
    help='a list of reads names')
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

    # read read names
    reads = read_readNames(args.reads)
    print("done reading read names:",args.reads)

    # read bam files
    bam1 = pysam.AlignmentFile(args.bam1, 'rb')
    bam2 = pysam.AlignmentFile(args.bam2, 'rb')

    # index by read name
    print("indexing file by read name:",args.bam1)
    name_indexed1 = index_by_readName(bam1)
    print("indexing file by read name:",args.bam2)
    name_indexed2 = index_by_readName(bam2)

    # compare reads
    comp_map(bam1, bam2, name_indexed1, name_indexed2, reads, args.base, args.suffix1, args.suffix2, args.AS_tag, args.NM_tag)
    print("done")

# to read read names
def read_readNames(file):
    with open(file, "r") as f:
        reads = f.read().splitlines()
    return reads

# index bam by read name
def index_by_readName(bam):
    name_indexed = pysam.IndexedReads(bam)
    name_indexed.build()
    return name_indexed

# function to compare reads
def comp_map(bam1, bam2, name_indexed1, name_indexed2, reads, outbase, suffix1, suffix2, as_tag, nm_tag):
    # count alignments
    def count_aln(aln):
        c = 0
        for a in aln:
            c += 1
        return c

    # return score and mismatches
    def get_as_nm(x, as_tag, nm_tag):
        return x.get_tag(as_tag), x.get_tag(nm_tag)

    def write_aln(x, out, tag):
        x.tags += [('CM',tag)]
        out.write(x)

    def compare_and_write(x1, x2, as_tag, nm_tag, out1, out2):
        score1, mismatches1 = get_as_nm(x1, as_tag, nm_tag)
        score2, mismatches2 = get_as_nm(x2, as_tag, nm_tag)
        if score1 > score2:
            write_aln(x1, out1, 'M')
        elif score1 < score2:
            write_aln(x2, out2, 'M')
        else:
            if mismatches1 < mismatches2:
                write_aln(x1, out1, 'M')
            elif mismatches1 > mismatches2:
                write_aln(x2, out2, 'M')
            else:
                write_aln(x1, out1, 'A')
                write_aln(x2, out2, 'A')

    # copy header
    header1 = bam1.header.copy()
    header2 = bam2.header.copy()
    # output files
    out1 = pysam.Samfile(outbase+"_"+suffix1+".bam", 'wb', header=header1)
    out2 = pysam.Samfile(outbase+"_"+suffix2+".bam", 'wb', header=header2)
    #out_amb1 = pysam.Samfile(outbase+"_amb_"+suffix1+".bam", 'wb', header=header1)
    #out_amb2 = pysam.Samfile(outbase+"_amb_"+suffix2+".bam", 'wb', header=header2)
    counter = 0
    for name in reads:
        counter += 1
        if name[0] == "@":
            name = name[1:]
        name = name.split(" ")[0]
        try:
            x = name_indexed1.find(name)
        except KeyError:
            try:
                x = name_indexed2.find(name)
            except KeyError:
                pass
            else:
                if count_aln(x) == 1:
                    x = next(x)
                    if not any([x.is_unmapped, x.is_secondary, x.is_supplementary]):
                        write_aln(x, out2, 'M')
        else:
            try:
                x = name_indexed2.find(name)
            except KeyError:
                if count_aln(x) == 1:
                    x = next(x)
                    if not any([x.is_unmapped, x.is_secondary, x.is_supplementary]):
                        write_aln(x, out1, 'M')
            else:
                # count alignments for each read
                unique_matches1 = count_aln(name_indexed1.find(name))
                unique_matches2 = count_aln(name_indexed2.find(name))
                # check if alignments are unique
                if unique_matches1 == 1 and unique_matches2 == 1:
                    x1 = next(name_indexed1.find(name))
                    x2 = next(name_indexed2.find(name))
                    # test if reads are mapped
                    if not any([x1.is_unmapped, x1.is_secondary, x1.is_supplementary]) and not any([x2.is_unmapped, x2.is_secondary, x2.is_supplementary]):
                        compare_and_write(x1, x2, as_tag, nm_tag, out1, out2)
                    elif not any([x1.is_unmapped, x1.is_secondary, x1.is_supplementary]) and any([x2.is_unmapped, x2.is_secondary, x2.is_supplementary]):
                        write_aln(x1, out1, 'M')
                    elif any([x1.is_unmapped, x1.is_secondary, x1.is_supplementary]) and not any([x2.is_unmapped, x2.is_secondary, x2.is_supplementary]):
                        write_aln(x2, out2, 'M')
                elif unique_matches1 == 1 and unique_matches2 > 1:
                    x1 = next(name_indexed1.find(name))
                    if not any([x1.is_unmapped, x1.is_secondary, x1.is_supplementary]):
                        write_aln(x1, out1, 'M')
                elif unique_matches1 > 1 and unique_matches2 == 1:
                    x2 = next(name_indexed2.find(name))
                    if not any([x2.is_unmapped, x2.is_secondary, x2.is_supplementary]):
                        write_aln(x2, out2, 'M')
        if counter % 100000 == 0:
            print("parsed "+str(counter)+" records")
    out1.close()
    out2.close()

if __name__ == "__main__":
    main()
