#!/bin/bash
# This bash script can extract protein fasta file based on screened busco file
# Then generate per-BUSCO fasta with one sequence per species

TSV_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/filtered_tsv
FAA_DIR=~/Binp29_exercise/MalariaProject/4_findortho
FILTERED_FAA_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/busco_sequence
OUT_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/per_busco_faa

mkdir -p "$FILTERED_FAA_DIR"
mkdir -p "$OUT_DIR"

species=(pf pb pc pk pv py tg ht)

# step1: extract filtered faa
for file in "$TSV_DIR"/*.tsv; do
	base=$(basename "$file" _full_table_filtered.tsv)

	awk -F'\t' '{print $3}' $file > "${TSV_DIR}/${base}_ids.txt"
	seqkit grep -f "${TSV_DIR}/${base}_ids.txt" "${FAA_DIR}/${base}.faa" -o "${FILTERED_FAA_DIR}/${base}_filtered.faa"
	echo "Filtered fasta for $base done"
done

# step2: make per-busco fasta
awk -F'\t' '{print $1}' "$TSV_DIR/ht_full_table_filtered.tsv" | sort | uniq > "$TSV_DIR/busco_list.txt"

while read busco; do
    out_file="$OUT_DIR/${busco}.faa"
    > "$out_file"

    for sp in "${species[@]}"; do
	    seq_id=$(awk -F'\t' -v b="$busco" '$1==b {print $3}' "$TSV_DIR/${sp}_full_table_filtered.tsv")
	    
	    if [ -z "$seq_id" ]; then
            echo "Warning: BUSCO $busco missing in $sp"
            continue
        fi


	seq=$(seqkit grep -r -p "^${seq_id}$" "$FILTERED_FAA_DIR/${sp}_filtered.faa" | sed "s/^>.*/>${sp}/")	
	echo "$seq" >> "$out_file"

	done
done < "$TSV_DIR/busco_list.txt"
