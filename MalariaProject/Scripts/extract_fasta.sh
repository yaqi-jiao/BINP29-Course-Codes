#!/usr/bin/env bash

GENOME_DIR=~/Binp29_exercise/MalariaProject/1_raw/genome
GTF_DIR=~/Binp29_exercise/MalariaProject/1_raw/
SCRIPT_DIR=~/Binp29_exercise/MalariaProject/Scripts

for genome in "${GENOME_DIR}"/*.genome; do
	[ -f "$genome" ] || continue
		echo "processing $genome"
		base=$(basename ${genome} .genome)

		"${SCRIPT_DIR}"/gffParse.pl -c -p -i "$genome" -g "${GTF_DIR}/${base}.gtf" -b "${base}"
	done
