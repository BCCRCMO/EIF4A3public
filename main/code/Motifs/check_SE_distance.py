import sys
from argparse import ArgumentParser, FileType


def outputDistances(allmotifs, samplesDir, hitsDir):

    motifs = allmotifs.split(':')

    samples_pre = ['increasing', 'decreasing']
    samples_post = '_p_values'
    hits_post = '_outputs'
    inputs = [('exon_up', 'j_1', 'e'),    # (region, junction_number, distance from start or end?)
              ('up_5prime', 'j_1', 's'),
              ('up_3prime', 'j_2', 'e'),
              ('exon_se', 'j_2', 's'),
              ('exon_se', 'j_3', 'e'),
              ('dn_5prime', 'j_3', 's'),
              ('dn_3prime', 'j_4', 'e'),
              ('exon_dn', 'j_4', 's')
              ]

    for x in range(len(inputs)):
        #     for x in range(1):
        for pre in samples_pre:
            sampleFile = samplesDir + '/' + pre + '_' + inputs[x][0] + samples_post   # positive/negative samples
            hitsFile = hitsDir + '/merged_' + inputs[x][0] + hits_post   # motif hits and their locations
            motifhits = {}

            with open(hitsFile, 'r') as f:
                content = f.readlines()
                for cnt in content:
                    cols = cnt.rstrip().split('\t')
                    id = cols[0].split('_')[0]

                    if cols[1] in motifs:
                        start = int(cols[0].split('_')[3])
                        end = int(cols[0].split('_')[4])
                        dist = 1000000

                        if inputs[x][2] == 'e':
                            dist = -(end - start - int(cols[2]))
                        else:
                            dist = int(cols[2])

                        if abs(dist) <= 300:
                            if id in motifhits:
                                motifhits[id].append(dist)
                            else:
                                motifhits[id] = [dist]

            for key in motifhits:
                motifhits[key] = list(set(motifhits[key]))


            lens = {}
            with open(sampleFile, 'r') as f:
                content = f.readlines()
                lens['positive'] = len([1 for j in content if j.rstrip().split(' ')[0] == 'positive'])
                lens['negative'] = len([1 for j in content if j.rstrip().split(' ')[0] == 'negative'])

            with open(sampleFile, 'r') as f:
                content = f.readlines()
                for cnt in content:
                    cols = cnt.rstrip().split(' ')
                    if cols[1] in motifhits:
                        for m in motifhits[cols[1]]:
                            print inputs[x][1], pre, cols[0], cols[1], m, float(1)/lens[cols[0]]


if __name__ == '__main__':

    parser = ArgumentParser(description=
                            '''Plotting distribution of motif hits for a set of given motifs''')
    parser.add_argument('--hitsdir', type=str,
                        help='Motif hits directory')
    parser.add_argument('--samplesDir', type=str,
                        help='Directory of samples')
    parser.add_argument('--allmotifs', type=str,
                        help='Selected motifs, separated by :')

    args = parser.parse_args()
    outputDistances(args.allmotifs, args.samplesDir, args.hitsdir)
