#!/usr/bin/env bash

FAA_DIR=~/Binp29_exercise/MalariaProject/4_findortho

for faa in "${FAA_DIR}"/*.faa; do
        [ -f "$faa" ] || continue
                echo "processing $faa"
                base=$(basename ${faa} .faa)

                busco -i "${FAA_DIR}/${base}.faa" -o "${base}" -m prot -l apicomplexa
        done
