import sys
import pandas as pd
from Bio import SeqIO
from Bio.motifs.matrix import PositionWeightMatrix
from Bio.Alphabet import IUPAC

notEmptyFiles = "motif_file"
notEmpty = []

def calc_bg_freq(seqrecs):
    freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    nt = 0
    for seqrec in seqrecs:
        seq = seqrec.seq.upper()
        nt += len(seq)
        for b in 'ACGT':
            freq[b] += float(seq.count(b))
    for b in 'ACGT':
        freq[b] /= nt
    return freq


def read_pwm_data(pwm_filename):
    pwm = pd.read_table(pwm_filename)
    pwm = {'A': pwm['A'].tolist(),
           'C': pwm['C'].tolist(),
           'G': pwm['G'].tolist(),
           'T': pwm['U'].tolist()}
    return pwm


def get_pssms(pwm_info, pwm_dir, bg_freq):


    with open(notEmptyFiles, 'r') as f:
        inputLines = f.readlines()
        for x in xrange(len(inputLines)):
            notEmpty.append(inputLines[x].rstrip())

    #pwm_info = pwm_info[pwm_info['MSource_Author'] == 'Ray'].copy()            
    pwm_info = pwm_info[pwm_info['Motif_ID'].isin(notEmpty)].copy()
    pwm_info.sort(columns='Motif_ID', inplace=True)

    pssms = {}

    motif_ids = pwm_info['Motif_ID'].drop_duplicates()
    for motif_id in motif_ids:
        pwm_data = read_pwm_data(pwm_dir + '/' + motif_id + '.txt')
        pwm = PositionWeightMatrix(IUPAC.unambiguous_dna, pwm_data)
        pssm = pwm.log_odds(bg_freq)
        pssms[motif_id] = pssm, pssm.max * 0.8

    print pssms
    return pssms


def get_pssm_matches(pssm, seq, thr):
    matches = []
    for pos, score in pssm.search(seq, threshold=thr, both=False):
        matches.append([pos, score])
    if len(matches) > 0:
        return pd.DataFrame(matches, columns=['seq_pos', 'score'])
    else:
        return pd.DataFrame(columns=['seq_pos', 'score'])


def get_all_pssm_matches(seq, pssms):
    all_pssm_matches = []
    for motif_id in pssms.iterkeys():
        pssm, thr = pssms[motif_id]
        pssm_matches = get_pssm_matches(pssm, seq, thr)
        pssm_matches['motif_id'] = motif_id
        all_pssm_matches.append(pssm_matches)
    flag = False
    for df in all_pssm_matches:
        if not df.empty:
            flag = True
    if flag:
        return pd.concat(all_pssm_matches)
    else:
        return []


def write_seq_num_progress(n_seqs):
    out_str = 'analyzing sequence {}'.format(n_seqs)
    sys.stdout.write('\r')
    sys.stdout.write(out_str)
    sys.stdout.flush()
    return


def get_motif_matches(seq_file, pwm_info_file, pwm_dir):
    seqrecs = SeqIO.parse(seq_file, 'fasta', alphabet=IUPAC.unambiguous_dna)
    bg_freq = calc_bg_freq(seqrecs)

    bg_str = 'background freq: A:{:.2f}, C:{:.2f}, G:{:.2f}, T:{:.2f}'
    print bg_str.format(bg_freq['A'], bg_freq['C'], bg_freq['G'], bg_freq['T'])

    pwm_info = pd.read_table(pwm_info_file)
    pssms = get_pssms(pwm_info, pwm_dir, bg_freq)
    print 'motif count: {}'.format(len(pssms))

    matches = []
    seq_ctr = 1
    seqrecs = SeqIO.parse(seq_file, 'fasta', alphabet=IUPAC.unambiguous_dna)
        
    for seqrec in seqrecs:

        write_seq_num_progress(seq_ctr)

        seq_ctr += 1

        seq = seqrec.seq.upper()
        
        event_name = seqrec.description

        all_motif_matches = get_all_pssm_matches(seq, pssms)

        if len(all_motif_matches) > 0:
            all_motif_matches['event_name'] = event_name
            matches.append(all_motif_matches)
    print '\n'

    cols = ['event_name', 'motif_id', 'seq_pos', 'score']
    matches = pd.concat(matches)[cols]
    matches['seq_pos'] = matches['seq_pos'].astype(int)
    return matches


if __name__ == '__main__':
    from argparse import ArgumentParser

    desc = 'find PWM motif matches'
    p = ArgumentParser(description=desc)

    p.add_argument('seqs', help='fasta file of seqs to search')
    p.add_argument('pwm_info', help='tsv of motif information')
    p.add_argument('pwm_dir', help='path to directory containing PWMs')
    p.add_argument('matches', help='output tsv of motif matches')

    args = p.parse_args()

    matches = get_motif_matches(args.seqs, args.pwm_info, args.pwm_dir)
    matches.to_csv(args.matches, sep='\t', index=False)
