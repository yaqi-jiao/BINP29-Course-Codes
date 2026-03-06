#!/bin/bash
# This script is using clustao to process faa files in `per_busco_faa`. It create a list for all busco ids, then run clustalo to generate alignment

FAA_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/per_busco_faa
ALIGN_DIR=~/Binp29_exercise/MalariaProject/5_align
TREE_DIR=~/Binp29_exercise/MalariaProject/6_tree

mkdir -p "$ALIGN_DIR"
mkdir -p "$TREE_DIR"

for faa in "${FAA_DIR}"/*; do
	
	base=$(basename "$faa" .faa)
	echo "Processing busco id $base"

	clustalo -i "$faa" -o "${ALIGN_DIR}"/"${base}"_aligned.faa --threads 10 -v
	raxmlHPC -s "${ALIGN_DIR}"/"${base}"_aligned.faa -n "${base}" -w "${TREE_DIR}" -o tg -m PROTGAMMABLOSUM62 -p 12345
done

