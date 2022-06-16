#!/bin/bash

for next in $1
do
	echo "$(cut -f2 $next)" | taxonkit lineage | taxonkit reformat -P | csvtk -H -t cut -f 1,3;
done > tmp_lca2
