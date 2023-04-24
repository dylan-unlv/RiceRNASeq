#!/bin/bash

COLORS=$(cat data/gene_modules.nofilter.txt | cut -f2 | sort | uniq)

for C in $COLORS
do
	grep $C data/gene_modules.nofilter.txt | cut -f1 > data/${C}_genes.nofilter.txt
done

