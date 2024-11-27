# Gut Microbial Mimics of Human Signaling Molecules

## Initial setup


1. Uniprot accession for human signalling proteins were retreived by filtering for 
 a. signal transduction (GO:0007165)
 b. signalling receptor activator activity (GO:0030546)
 c. extracellular space (GO:0005615)
 d. human taxanomy (9606)
 e. swissprot reviewed

2. Alphafold Structures for these accessions were retreived
3. Uniprot Accessions for  metagenomic proteins from the waldron dataset were retreived
4. Alphafold Structures for tehse accessions were retreived

```bash

curl "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28go%3A0007165%29+AND+%28go%3A0030546%29+AND+%28go%3A0005615%29+AND+%28taxonomy_id%3A9606%29+AND+%28reviewed%3Atrue%29%29" -o human_extracell_accs
cat human_extracell_accs | parallel "curl 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb' -o human_extracell_structures/{}.pdb"
grep ">" /workdir1/shared/20241121_all_gut_prot_unire90_repr.fasta | awk -F"_" '{print $2}' > waldron_accs
cat waldron_accs | parallel "curl 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb' -o waldron_structures/{}.pdb"

```



## Foldseek search and quality filtering

1. Folseek easy-search
2. The dataset was filtered for a foldseek probability score of 1, an alignment TM score of >= 0.5, and a query coverage of >= 80%
4. Taxanomic lineages for metegenomic proteins were retreievd using uniprots ID mapping tool
5. Non-eukaryotic accessions were removed
6. The dataframe was filtered for  non-eukaryotic sequences





```bash

nohup  ~/foldseek/bin/foldseek easy-search human_extracell_structures waldron_structures extracell_hits.m8 tmpfolder --threads 50 --format-output "query,target,pident,alntmscore,prob,evalue,bits,alnlen,qstart,qend,tstart,tend,qseq,tseq,qcov,tcov" &
awk -F"\t" '$5 == 1 && $4 >= 0.5 && $15 >= 0.8 {print $0}' extracell_hits.m8 > extracell_hits_prob1_TM50_qcov80.m8
curl "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/7ce976bbf7a6e7aa0a8694b71de7c98d73d5ab59?fields=accession%2Clineage&format=tsv" -o baccs.u.tax.tsv
grep "Archaea\|Bacteria"   baccs.u.tax.tsv  >  baccs.u.tax.bact.arch.tsv
awk -F"\t" 'FNR==NR {arr[$1] ; next} $2 in arr {print $0}'  baccs.u.tax.bact.arch.tsv   extracell_hits_prob1_TM50_qcov80.m8 > extracell_hits_prob1_TM50_qcov80_bact_arch.m8

```

## Signal Peptide Searching


1. The sequences of the metagenomic proteins where put into a fasta file
2. DeepTMHMM was run on those sequences
3. The signal peptide containing accessions were extracted
4. The alignmnet data was subsetted for thoe alignments were the metagenomic protein had  a signal peptide.

```bash
awk -F"\t" '{if (!seen[$2]++) {print ">" $2 "\n" $14}}' extracell_hits_prob1_TM50_qcov80_bact_arch.m8 > extracell_hits_prob1_TM50_qcov80_bact_arch.fa
docker1 run -v /workdir/djl294/:/openprotein/data/ -w /openprotein -e LC_ALL=C.UTF-8 -e LANG=C.UTF-8 --rm 7eb681f7f663 python3 predict.py --fasta data/extracell_hits_prob1_TM50_qcov80_bact_arch.fa
grep "| SP"  deeptmhmm_results/predicted_topologies.3line | awk -F" |" '{print $1}'  | awk -F">" '{print $2}' > sp_accs
awk -F"\t" 'FNR==NR {arr[$1]; next} $2 in arr {print $0}' sp_accs  extracell_hits_prob1_TM50_qcov80_bact_arch.m8 > extracell_hits_prob1_TM50_qcov80_bact_arch.sp.m8
```


