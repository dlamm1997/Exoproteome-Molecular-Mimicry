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
2. The dataset was filtered for a foldseek probability score of 1 and alignment TM score of >= 0.5
3. Taxanomic lineages for metegenomic proteins were retreievd using uniprots ID mapping tool
4. Non-eukaryotic accessions were removed
5. The dataframe was filtered for  non-eukaryotic sequences





```bash

nohup  ~/foldseek/bin/foldseek easy-search human_extracell_structures waldron_structures extracell_hits.m8 tmpfolder --threads 50 --format-output "query,target,pident,alntmscore,prob,evalue,bits,alnlen,qstart,qend,tstart,tend,qseq,tseq,qcov,tcov" &
awk -F"\t" '$5 == 1 && $4 >= 0.5 {print $0}' extracell_hits.m8 > extracell_hits_prob1_TM50.m8
curl "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/230df76ad8570db60e45823331400d08d084fa4b?fields=accession%2Clineage&format=tsv" -o baccs.u.tax.tsv
grep -v "Eukaryota" baccs.u.tax.tsv > baccs.u.tax.not_euk
awk -F"\t" 'FNR==NR {arr[$1] ; next} $2 in arr {print $0}'  baccs.u.tax.not_euk extracell_hits_prob1_TM50.m8  > extracell_hits_prob1_TM50.not_euk.m8

```

## Signal Peptide Searching


1. The sequences of the metagenomic proteins where put into a fasta file
2. signalp5 was run for both gram+ and gram- 
3. The signalp outputs were sorted for those that have some sort of signal peptide
4. After combining those files, and taking the unique list of proteins, that alignment dataframe was further filtered for only those with postive signal peptide

```bash
awk -F"\t" '{if (!seen[$2]++) {print ">" $2 "\n" $14}}' extracell_hits_prob1_TM50.m8 > extracell_hits_prob1_TM50.bact.fa
nohup /programs/signalp-5.0/bin/signalp -org gram- -fasta extracell_hits_prob1_TM50.not_euk.bact.fa -gff3  &
awk -F"\t" 'NR > 2 && !($2 ==  "OTHER") {print $1}' gram_pos/extracell_hits_prob1_TM50.not_euk.bact_summary.signalp5 > gram_pos_sig ; awk -F"\t" 'NR > 2 && !($2 ==  "OTHER") {print $1}' gram_neg/extracell_hits_prob1_TM50.not_euk.bact_summary.signalp5 > gram_neg_sig 
awk -F"\t" 'FNR==NR {arr[$1] ; next} $2 in arr {print $0}' all_sig.u extracell_hits_prob1_TM50.not_euk.m8 > extracell_hits_prob1_TM50.not_euk.sp5.m8
```


# Additional Filters

1. We filtered for query coverage > 80%

```bash
awk -F"\t" '$15 > 0.8 {print $0}' extracell_hits_prob1_TM50.not_euk.sp5.m8 > extracell_hits_prob1_TM50.not_euk.sp5.qcov80.m8
```
