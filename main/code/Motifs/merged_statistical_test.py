import sys
import random
from sets import Set
from scipy.stats import mannwhitneyu
from scipy.stats import rankdata
from scipy.stats import kstest
from scipy.stats import ks_2samp

import numpy as np
from scipy import stats

motifCounts = {}
all_motifs = []
p_values = {}
num_wins = {}
positives, backgrounds = {}, {}
event_len = {}
positive_nonZeros = {}
negative_nonZeros = {}
positive_hit_ratio = {}
negative_hit_ratio = {}
parentMotif = {}
counted = {}                                   
normalizedCounts = {}


def load_info(info_file):
    with open(info_file, 'r') as f:
        content = f.readlines()
        for line in content:
            cols = line.rstrip().split(' ')
            event_len[cols[0]] = int(cols[4]) - int(cols[3]) + 1
    
def weighted_random(choices, num):

    out = []
    for k in xrange(num):
#         print k
        total = sum(w for c, w in choices)
        r = random.uniform(0, total)
        curr = 0
        # print len(choices)
        for i in xrange(len(choices)):
            if curr + choices[i][1] >= r:
                out.append(choices[i][0])
                break
            curr += float(choices[i][1])
            if i == len(choices)-1:
                raise Exception("Error in random number generation!")

    return out



def count_hits(motifFile, mergedMotifsFile):

    mergedMotifsInput = mergedMotifsFile.readlines()
    for i in xrange(len(mergedMotifsInput)):
        cols = mergedMotifsInput[i].rstrip().split(' ')
        for j in xrange(len(cols)):
            parentMotif[cols[j]] = '_'.join(cols)
    
    motifInput = motifFile.readlines()
    global all_motifs
    
    for i in xrange(len(motifInput)):

        cols = motifInput[i].split('\t')
        eventId, motifId, loc = cols[0].split('_')[0], parentMotif[cols[1]], int(cols[2])
        all_motifs.append(parentMotif[cols[1]])

        s = (eventId, motifId, loc)
        t = (eventId, motifId)

        if (not (s in counted)) and eventId in event_len:
            counted[s] = 1
            if t in motifCounts:
                motifCounts[t] = motifCounts[t] + 1
            else:
                motifCounts[t] = 1

    for key in motifCounts:

        #start, end = int(key[0].split('_')[3]), int(key[0].split('_')[4])
        motifCounts[key] = float(motifCounts[key]) / event_len[key[0]]
        #normalizedCounts[(key[0].split('_')[0], key[1])] = motifCounts[key] / float(end-start+1)

    #     for key in motifCounts:
    #         print key, motifCounts[key]
    all_motifs = list(Set(all_motifs))


def statistical_test(positive_file, event_type):

    positives, negatives = [], []
    with open(positive_file, 'r') as f:
        content = f.readlines()
        for x in content:
            if x.rstrip() in event_len:
                positives.append(x.rstrip())

    for key in event_len:
        if not (key in positives):
            negatives.append(key)

    negative_intervals = []
    for i in xrange(len(negatives)):
        negative_intervals.append((event_len[negatives[i]], negatives[i]))
    negative_intervals.sort()    # sort negative samples based on their lengths

#     for i in range(len(negative_intervals)):
#         print negative_intervals[i]
        
    for m in all_motifs:
        p_values[m] = []
        num_wins[m] = []
        positive_nonZeros[m] = []
        negative_nonZeros[m] = []
        
    # SPLIT BACKGROUND BASED ON THEIR LENGTHS
    NUM = 30       # Number of groups to consider
    lows = []
    highs = []
    counts = [0]*NUM
    group = [0]*len(negative_intervals)
        
    curr = 0      # define low and high counts for each quantile
    for i in xrange(NUM):  # lows and highs define intervals' starts and ends
        high = int((len(negative_intervals)*(i+1))/NUM)-1
        lows.append(curr)
        highs.append(high)
        for j in xrange(curr, high+1):
            group[j] = i
        curr = high+1

    
    for i in xrange(len(positives)):
        curLen = event_len[positives[i]] 

        if curLen < negative_intervals[0][0]:             # In first group
            counts[0] += 1
        elif curLen > negative_intervals[-1][0]:          # In last group
            counts[-1] += 1
        else:
            candidates = []
            for j in xrange(NUM):
                if (curLen >= negative_intervals[lows[j]][0] and curLen <= negative_intervals[highs[j]][0]):    # Find all intervals that a positive sample belongs to
                    candidates.append(j)
            for x in candidates:
                counts[x] += (float(1) / len(candidates))                   # Split the sample between them

    # counts contains the number of samples to be drawn from each interval for negative samples
    #negatives = []
    FACTOR = 10            # Ratio of negative to positive samples
    num_samples = int(FACTOR*len(positives))
    #     print num_samples
    choices = []
    for i in xrange(len(negative_intervals)):            # Assign weight in each bin based on number of positive samples in that bin
        if counts[group[i]] > 0:
            choices.append((negative_intervals[i][1], counts[group[i]]))

    #     negative_samples = negatives   # all intronic regions have the same length (e.g. 300)
    #     negative_samples = np.random.choice(negatives, num_samples, replace=False)
    #     if event_type == 'exon':
    negative_samples = weighted_random(choices, num_samples)   # Pick samples randomly
    positive_samples = positives
    
    for i in xrange(len(positive_samples)):
        print 'positive', positive_samples[i]
    for i in xrange(len(negative_samples)):
        print 'negative', negative_samples[i]

    #     for i in range(len(negative_samples)):
    #         print negative_samples[i]

    for motif in all_motifs:

        pValues = []
        nValues = []
        pos_lens = []
        neg_lens = []
        pos_counts = []
        neg_counts = []
        
        for i in xrange(len(positive_samples)): # Calculate positive hit ratio for each motif
            pValues.append(0)
            pos_counts.append(0)
            if (positive_samples[i], motif) in motifCounts:
                pValues[-1] = motifCounts[(positive_samples[i], motif)]
                pos_counts[-1] = motifCounts[(positive_samples[i], motif)] * event_len[positive_samples[i]]
            pos_lens.append(event_len[positive_samples[i]])
        
        for i in xrange(len(negative_samples)): # Calculate negative hit ratio for each motif
            nValues.append(0)
            neg_counts.append(0)
            if (negative_samples[i], motif) in motifCounts:
                nValues[-1] = motifCounts[(negative_samples[i], motif)]
                neg_counts[-1] = motifCounts[(negative_samples[i], motif)] * event_len[negative_samples[i]]
            neg_lens.append(event_len[negative_samples[i]])


        len1 = sum(x for x in pos_lens) / len(pos_lens)
        len2 = sum(x for x in neg_lens) / len(neg_lens)

        #   print len1, len2
        
        if event_type == 'exon' and (float(len1)/float(len2) > 1.15 or 
               float(len1)/float(len2) < 0.85):
            raise Exception("Samples did not split well!")

        #         if (float(len1)/float(len2) > 1.15 or 
        #             float(len1)/float(len2) < 0.85):
        #             raise Exception("Samples did not split well!")


        ranks = rankdata(pValues + nValues)     # This part is not used here
        ranks_x = ranks[0:len(pValues)]
        ranks_y = ranks[len(pValues):]

        wins = 0
        for k in xrange(10000):
            x = random.sample(xrange(len(pValues)), 1)[0]
            y = random.sample(xrange(len(nValues)), 1)[0]
            if pValues[x] > nValues[y]:
                wins += 1
            elif pValues[x] < nValues[y]:
                wins -= 1

        pValues_positive_ratio = float(sum(i > 0 for i in pValues)) / len(pValues)
        nValues_positive_ratio = float(sum(i > 0 for i in nValues)) / len(nValues)

        # print motif, float(sum(i > 0 for i in pValues)), len(pValues)
        # print motif, float(sum(i > 0 for i in nValues)), len(nValues)

        positive_nonZeros[motif].append(pValues_positive_ratio)
        negative_nonZeros[motif].append(nValues_positive_ratio)

        positive_hit_ratio[motif] = sum(pos_counts) / sum(pos_lens)
        negative_hit_ratio[motif] = sum(neg_counts) / sum(neg_lens)

        if(positive_nonZeros[motif][-1] == 0 and \
           negative_nonZeros[motif][-1] == 0):
            p_values[motif].append(0.5)
        else:
            p_values[motif].append(mannwhitneyu(pValues, nValues)[1])
        num_wins[motif].append(wins)


    for motif in all_motifs:
        if np.mean(np.array(num_wins[motif])) > 0:
            print motif, np.mean(np.array(p_values[motif])),
        else:
            print motif, -np.mean(np.array(p_values[motif])),
        print np.mean(np.array(positive_nonZeros[motif])), \
              np.mean(np.array(negative_nonZeros[motif])), \
              positive_hit_ratio[motif], \
              negative_hit_ratio[motif], \
              positive_hit_ratio[motif] / (negative_hit_ratio[motif]+0.0000000001)


if __name__ == '__main__':

    from argparse import ArgumentParser, FileType

    parser = ArgumentParser(description=
                            '''Performs statistical tests to find enriched motifs''')
    parser.add_argument('--motifs', type=FileType('r'),
                        help='Motif hit file')
    parser.add_argument('--mergedMotifs', type=FileType('r'),
                        help='Guide on how to merge motifs')
    parser.add_argument('--type', type = str,
                        help = 'Type of region (intorn | exon)')
    parser.add_argument('--positive', type = str,
                        help = 'Positive ids')
    parser.add_argument('--events_info', type = str,
                        help = 'start, end, strand information for events')
    args = parser.parse_args()

    load_info(args.events_info)
    count_hits(args.motifs, args.mergedMotifs)
    statistical_test(args.positive, args.type)
