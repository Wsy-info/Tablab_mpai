import scrublet as scr
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

### sample1
infile = "sample1_count.txt"
outfile = "sample1_Scrublet_result.txt"
finallist = []
with open(infile, 'r') as f:
    header = next(f)
    cell_barcodes = header.rstrip().split('\t')
    for line in f:
        tmpline = line.rstrip().split('\t')[1: ]
        tmplist = [float(s) for s in tmpline]
        finallist.append(tmplist)

finalarray = np.array(finallist)
count_matrix = np.transpose(finalarray)

scrub = scr.Scrublet(count_matrix, expected_doublet_rate = 0.08)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.savefig('sample1/doublet_histogram.pdf', dpi = 300)
predicted_doublets_final = scrub.call_doublets(threshold = 0.25)

with open(outfile, 'w') as f:
    f.write('\t'.join(['CB', 'Scrublet', 'Scrublet_Score']) + '\n')
    for i in range(len(doublet_scores)):
        if predicted_doublets_final[i] == 0:
            result = 'Singlet'
        else:
            result = 'Doublet'
            f.write('\t'.join([cell_barcodes[i], result, str(doublet_scores[i])]) + '\n')
