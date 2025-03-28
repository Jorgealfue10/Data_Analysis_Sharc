#!/bin/bash

file=QM.out
fileout=QM_test.out
for dir in ICOND*
do
	if [ -d $dir ] 
	then
		echo $dir
		cd $dir
		cp ../dmbydy.py .
		cp ../lib_QMout.py .
		if [ -f $file ] 
		then
			mkdir bckup
			cp bckup/QM.out $file 
			cp $file bckup/.
			python3 dmbydy.py
			cp $fileout $file
			cp $fileout bckup/.
		else
			echo "File not found. " $dir
		fi
		cd -
	fi
done
