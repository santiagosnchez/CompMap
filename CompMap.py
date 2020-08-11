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
python CompMap.py -1 first.bam -2 second.bam -r read_name_list -b my_base_name

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
    compare_reads(bam1, bam2, name_indexed1, name_indexed2, reads, args.base)
    print("done")

# to read read names
def read_readNames(file):
    with open(file, "r") as f:
        reads = f.readlines()
    return reads

# index bam by read name
def index_by_readName(bam):
    name_indexed = pysam.IndexedReads(bam)
    name_indexed.build()
    return name_indexed

# function to compare reads
def compare_reads(bam1, bam2, name_indexed1, name_indexed2, reads, outbase):
    # copy header
    header1 = bam1.header.copy()
    header2 = bam2.header.copy()
    # output files
    out1 = pysam.Samfile(outbase+"_matched1.bam", 'wb', header=header1)
    out2 = pysam.Samfile(outbase+"_matched2.bam", 'wb', header=header2)
    out_amb1 = pysam.Samfile(outbase+"_ambiguous1.bam", 'wb', header=header1)
    out_amb2 = pysam.Samfile(outbase+"_ambiguous2.bam", 'wb', header=header2)
    counter = 0
    for name in reads:
        counter += 1
        name = name.rstrip()
        if name[0] == "@":
            name = name[1:]
        name = name.split(" ")[0]
        try:
            name_indexed1.find(name)
        except KeyError:
            try:
                name_indexed2.find(name)
            except KeyError:
                pass
            else:
                for x in name_indexed2.find(name):
                    out2.write(x)
        else:
            try:
                name_indexed2.find(name)
            except KeyError:
                for x in name_indexed1.find(name):
                    out1.write(x)
            else:
                score1 = 0
                score2 = 0
                mismatches1 = 0
                mismatches2 = 0
                for x in name_indexed1.find(name):
                    if not x.is_unmapped or not x.is_secondary or not x.is_supplemenrary:
                        score1 = x.get_tag(args.AS_tag)
                        mismatches1 = x.get_tag(args.NM_tag)
                for x in name_indexed2.find(name):
                    if not x.is_unmapped or not x.is_secondary or not x.is_supplemenrary:
                        score2 = x.get_tag(args.AS_tag)
                        mismatches2 = x.get_tag(args.NM_tag)
                if score1 > score2:
                    for x in name_indexed1.find(name):
                        out1.write(x)
                elif score2 > score1:
                    for x in name_indexed2.find(name):
                        out2.write(x)
                else:
                    if mismatches1 < mismatches2:
                        for x in name_indexed1.find(name):
                            out1.write(x)
                    if mismatches1 > mismatches2:
                        for x in name_indexed2.find(name):
                            out2.write(x)
                    else:
                        for x in name_indexed1.find(name):
                            out_amb1.write(x)
                        for x in name_indexed2.find(name):
                            out_amb2.write(x)
                    #it = name_indexed1.find(name)
                    #x = it.next()
                    #o.write("@"+name+"\n"+x.seq+"\n+\n"+x.qual+"\n")
        if counter % 100000 == 0:
            print("parsed",counter,"records")


if __name__ == "__main__":
    main()
