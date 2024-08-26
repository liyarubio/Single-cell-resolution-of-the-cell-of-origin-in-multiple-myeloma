#bash /media/liyaru/LYR16T/Blood/Code/2_Shell/8_SCENIC

i=CD19
pyscenic grn \
--num_workers 10 \
--output adj.${i}.tsv \
--method grnboost2 \
${i}.loom /home/liyaru/public_Data/SCENIC/hs_hgnc_tfs.txt

i=CD19
pyscenic ctx \
adj.${i}.tsv \
/home/liyaru/public_Data/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /home/liyaru/public_Data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ${i}.loom \
--mode "dask_multiprocessing" \
--output ${i}.reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
${i}.loom \
${i}.reg.csv \
--output ${i}.SCENIC.loom \
--num_workers 10
