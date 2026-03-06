#
# Description: this script is used to screen the busco results for full.tsv

FILE_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/tsv_set
OUT_DIR=~/Binp29_exercise/MalariaProject/4_findortho/BUSCO/filtered_tsv

mkdir -p "$OUT_DIR"

for file in "$FILE_DIR"/*; do
	base=$(basename "$file" .tsv)
	out="$OUT_DIR/${base}_filtered.tsv"
	
	awk -F '\t' '
	$2=="Complete" || $2=="Duplicated" {
	
		busco=$1
		score=$4
		
		if (!(busco in best) || score > best[busco]) {
		best[busco] = score
		line[busco] = $0
		}

	}
	
	END {
		for (b in line)
			print line[b]
	}
	
	' "$file" > "$out" | sort -k1,1 > "$out"
	echo "Processed $file → $out"
done
