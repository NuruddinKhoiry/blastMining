#!/bin/bash
for next in $1
do
	echo "$(cut -f2 $next)" | taxonkit lca  -s ","
done >> tmp_lca
