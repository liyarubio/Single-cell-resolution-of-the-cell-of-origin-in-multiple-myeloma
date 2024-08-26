cd /media/liyaru/LYR/MM2023/5_Result/17_scATAC
bedtools getfasta -fi /home/liyaru/public_Data/refdata-gex-GRCh38-2020-A/fasta/genome.fa -bed peak_FCRL5.bed > peak_FCRL5.bed.fa


cd ~/Software/Software/meme-5.1.1/src
./fimo

./fimo --oc /media/liyaru/LYR/MM2023/5_Result/17_scATAC/peak_FCRL5_fimo ../motif/JASPAR2022_CORE_non-redundant_pfms_meme.txt /media/liyaru/LYR/MM2023/5_Result/17_scATAC/peak_FCRL5.bed.fa
