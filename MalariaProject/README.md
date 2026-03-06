# Phylogenetic analysis of Malaria parasites

Research question: *Plasmodium falciparum* is more close to other mammalian parasites, or it's more close to a bird malaria parasite.

## Structure

```text
.
├── 1_raw: genome data and gtf files of *Plasmodium* species
│ 
├── 2_filter: first filter and removal
│  
├── 3_prediction: the final filtered Ht genome and gtf data
│ 
├── 4_findortho: find ortho proteins in 2 approaches
│   ├── orthoprotein
│   └── BUSCO: the final pre-busco faa files
│ 
├── 5_align: the result of alignment of each busco id
│ 
├── 6_tree: trees for each busco id, sample '*11330at5794' is provided
│ 
├── Script/
│
└── README.md
```

## workflow

## Step 1: Set-up directory

Create project directory `MalariaProject`, and create subdirectories.

Note: all process is based in this root directory

```bash
mkdir MalariaProject
cd MalariaProject/
mkdir 1_raw 2_filter 3_prediction 4_findortho 5_align 6_tree Scripts
```

## Step2: Collect data

Download data from server, including genome data of Haemoproteus_tartakovskyi, and gtf files of P species

```bash
cp /resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz ./1_raw/
cp /tmp/PutGenomeHere/* 1_raw/
```

## Step3:  Clean the data

Because the data derives both from the bird and the parasite, and only parasite data is required, the bird genome data is removed. There is a difference of GC content between bird genome(41%) and parasite genome(19%~42%).
The filter is set as follows:

    - GC content threshold:35%
    - scaffolds length < 3000 bp

As a result, sequence less than 3000 bp, and GC content larger than 35% will be removed using python script `RemoveScaffold.py`

```bash
./Scripts/RemoveScaffold.py ./1_raw/Haemoproteus_tartakovskyi.raw.genome 2_filter/filtered_H_t_genome
```
  
## Step4: Prediction

Use gmes.pl to generate gtf file from filtered genome data, then transfer gtf file into gff format.

```bash
/usr/local/bin/gmes_petap.pl --sequence 2_filter/filtered_Ht.genome --cores 10 --work_dir 3_prediction/ --ES --min_contig 5000
mv genemark.gtf Haemoproteus.gtf
cat Haemoproteus.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gtf
```

## Step5: Refilter

use `gff.Parse.pl`

```bash
../Scripts/gffParse.pl -c -p -i ../2_filter/filtered_Ht.genome -g Ht2.gtf
mkdir re-filter
cd re-filter
blastp -query ../Ht_protein_fasta -db SwissProt -out Ht_blastp.out -outfmt 6 -num_threads 10
python ../../Scripts/datParser.py Ht_blastp.out ../Ht.faa ../../1_raw/taxonomy.dat ../../1_raw/uniprot_sprot.dat > scaffolds.txt
python ../../Scripts/filter_fasta.py scaffolds.txt ../Ht.fna filtered_Ht.fna  # WROGN!!!
# A really scary problem is that I started confusing genome and FNA right from this step. Here, "filter" means to filter genome data, not FNA!!!
python ../../Scripts/filter_fasta.py scaffolds.txt ../../2_filter/filtered_Ht.genome filtered_Ht.genome
```
Check the result, then do the prediction again, and re-extract the faa sequence.

```bash
gmes_petap.pl --sequence ../filtered_Ht.genome --cores 12 --work_dir . --ES --min_contig 5000
mv genemark.gtf ht.gtf

```

## Step6: Find orthologs

extract fasta sequence of another species, bastch processing using bash script `extract_fasta.sh`

```bash
../../Scripts/extract_fasta.sh
```

1. use protein ortho
```bash
sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' *.faa
nohup proteinortho6.pl {ht,pb,pc,pf,pk,pv,py,tg}.faa &
```

2. use busco

```bash
nohup ../../Scripts/run_busco.sh &
```

After running, only one-to-one orthologs are kept. bash script `` is used to generate screened fasta files for each busco, then extract faa file for each id

```bash
mkdir tsv_set
cd tsv_set
for i in ht pb pc pf pk pv py tg; do
  mv ../"$i"/run_apicomplexa_odb12/full_table.tsv ../"$i"/run_apicomplexa_odb12/"$i"_full_table.tsv
  cp ../"$i"/run_apicomplexa_odb12/"$i"_full_table.tsv .
done
../../../Scripts/screen_buscos.sh
../../../Scripts/extract_proteinseq.sh
```

## Step7: alignment

use clustao to aligh all faa sequences using bash Script `multi_align`,
then run raxml for all alignments

```bash
../Scripts/multi_align.sh
```

## Step8: 
