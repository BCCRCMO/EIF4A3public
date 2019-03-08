import sys
import re

strand = {}
mRNAs = {}
lastExon = {}
geneName = {}
PSI_for = {} # either for proximal or distal

def extract_proximal_distal(ALEgff, direction, events):

    with open(ALEgff, 'r') as f:
        inputLines = f.readlines()[1:]
        for i in xrange(len(inputLines)):
            cols = inputLines[i].split('\t')
            l, r = int(cols[3]), int(cols[4])
            s = cols[6]
            if cols[2] == 'gene':
                gene_id = re.search('(?<=;gid=)\S+(?=\n)', cols[8]).group()
                strand[gene_id] = cols[6]
            elif cols[2] == 'mRNA':
                mRNA_id = re.search('(?<=ID=)\S+(?=;Parent)', cols[8]).group()
                gene_id = re.search('(?<=;gid=)\S+(?=\n)', cols[8]).group()
                strand[gene_id] = cols[6]
                mRNAs.setdefault(gene_id, []).append(mRNA_id)
                geneName[mRNA_id] = gene_id
            elif cols[2] == 'exon':
                mRNA_id = re.search('(?<=Parent=)\S+(?=;Name)', cols[8]).group()
                if not mRNA_id in lastExon:
                    lastExon[mRNA_id] = (l, r)
                else:
                    if (s == '+' and r > lastExon[mRNA_id][1]) or \
                       (s == '-' and l < lastExon[mRNA_id][0]):
                        lastExon[mRNA_id] = (l, r)

    for g in mRNAs:
        if len(mRNAs[g]) == 2:
            if (strand[g] == '+' and lastExon[mRNAs[g][0]][1] < lastExon[mRNAs[g][1]][1]) or \
               (strand[g] == '-' and lastExon[mRNAs[g][0]][0] > lastExon[mRNAs[g][1]][0]):
                PSI_for[g] = 'p'
            else:
                PSI_for[g] = 'd'


    with open(events, 'r') as f:
        inputLines = f.readlines()
        for i in xrange(len(inputLines)):
            id = inputLines[i].rstrip()
            if (direction == '+' and PSI_for[id] == 'p') or \
               (direction == '-' and PSI_for[id] == 'd'):
                print id, 'p'
            else:
                print id, 'd'


if __name__ == '__main__':
    from argparse import ArgumentParser, FileType

    parser = ArgumentParser(description='Find whether proximal or distal is overexpressed')
    parser.add_argument('--misoGff', type=str,
                        help='Miso gff file')
    parser.add_argument('--direction', type=str,
                        help='+ for monotonically increasing and - for decreasing')
    parser.add_argument('--events', type=str,
                        help='List of identified events')
    args = parser.parse_args()

    extract_proximal_distal(args.misoGff, args.direction, args.events)
