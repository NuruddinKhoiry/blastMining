#!/bin/bash
for next in $1
do
	echo "$(cut -f2 $next)" | taxonkit lca  -s "," -j $2;
done >> ./$3/TMPDIR/tmp_lca