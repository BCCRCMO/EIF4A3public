#!/Users/alborzmazloomian/anaconda/bin/python
import sys
import os

gene_sets = {}
all_genes = {}

with open('c5.bp.v5.2.symbols.gmt.txt', 'r') as f:
    content = f.readlines()
    for i in range(len(content)):
        cols = content[i].rstrip().split('\t')
        gene_sets[cols[0]] = {}
        for j in range(2, len(cols)):
            gene_sets[cols[0]][cols[j]] = 1
            all_genes[cols[j]] = 1

fout = open(sys.argv[1]+'.pathwaycounts', 'w')
total = len(all_genes.keys())
with open(sys.argv[1], 'r') as f:

    gene = [x.rstrip() for x in f.readlines()]
    for key in gene_sets:
        w_sample = 0
        t_sample = 0
        for i in range(len(gene)):
            if gene[i] in all_genes:
                t_sample += 1
            if gene[i] in gene_sets[key]:
                w_sample += 1
        fout.write(key+' '+str(w_sample)+' '+str(t_sample)+' '+str(len(gene_sets[key].keys()))+' '+str(total)+'\n')
        #print key, w_sample, t_sample, len(gene_sets[key].keys()), total
fout.close()

inputName = sys.argv[1]+'.pathwaycounts'
cmd = "./calc_fdr.R "+inputName+" > "+inputName+".stats"
os.system(cmd)
