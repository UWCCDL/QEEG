#!/bin/bash
# change csv format to tab format
delet_old=0
while getopts "d" arg
do
	case $arg in
		d) delet_old=1;;
		?) echo "Unknown argument"
		   exit 1;;
	esac
done


for file in $@;do
	file_basename=`basename $file .csv`
	absolute_path=`readlink -f $file`
	output_dir=`dirname absolute_path`
	cat $file | tr "," "\t" > $output_dir/$file_basename.txt
	echo Output:$output_dir/$file_basename.txt
	if [ $delet_old -eq 1 ];then
		rm -f $file
		echo $file has been deleted.
	fi
done
