for c in `ls | grep TR` 
do
	/catt/bwa-0.7.17/bwa index $c %c
done
