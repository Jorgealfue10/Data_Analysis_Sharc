#!/bin/bash

for i in ICOND*
do
	sed 's/SHARC_MOLPRO.py/parsing_MRCI.py/g' $i/run.sh > kk
	mv kk $i/run.sh
done
