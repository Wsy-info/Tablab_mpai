### Figure 4
### A

### pyscenic grn
pyscenic grn --num_workers 40 \
    --method grnboost2 \
    --seed 16 \
    --output brain.harmony.group.adj_new.tsv \
    brain.harmony_group_singlet_exp_new.loom \
    Reference/index_genome/cisTarget_databases/allTFs_mm.txt


### pyscenic ctx
pyscenic ctx --num_workers 40 \
    --output brain.harmony.group.regulons_new.csv \
    --expression_mtx_fname brain.harmony_group_singlet_exp_new.loom \
    --mask_dropouts \
    --mode "dask_multiprocessing" \
    --annotations_fname Reference/index_genome/cisTarget_databases/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
    --min_genes 10 \
    brain.harmony.group.adj_new.tsv \
    Reference/index_genome/cisTarget_databases/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

### pyscenic aucell
pyscenic aucell --num_workers 40 \
    --output brain.harmony.group_SCENIC_new.loom \
    brain.harmony_group_singlet_exp_new.loom \
    brain.harmony.group.regulons_new.csv


