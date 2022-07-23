#!/bin/bash

for next in $1
do
	echo "$(cut -f2 $next)" | taxonkit lineage -j $2 | taxonkit reformat -P -j $2 | csvtk -H -t cut -f 1,3 -j $2;
done > ./$3/TMPDIR/tmp_lca2