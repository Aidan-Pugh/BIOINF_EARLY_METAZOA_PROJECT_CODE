#!/bin/bash

# Run mafft on all files:
for i in *.fas
do
	mafft $i > $i.aln
done

## Trim the aligned sequences

# Load in trimAI
module load apps/trimAI/1.2

for i in *.fas.aln
do
	trimal -in $i -out $i.trimraw -automated1
done

exit
