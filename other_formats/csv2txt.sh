#!/bin/bash
# change csv format to tab format

for file in $@; do
	cat $file | tr "," "\t" 
done
