#!/bin/bash
#inputs:
#	$1	pdb dictionary
#	$2	exec file
#outputs:
#	{PDB DIR}/jd2_scores.sc
#	{PDB DIR}/{PDB}_0001.pdb
#	{PDB DIR}/score_jd2.log
#need:
#	NONE
pwd=`pwd`
touch $1/jd2_pdb.list
find $1 -name "*.pdb" |while read file
do
	dir=`dirname $file`
	cd $dir
	fpwd=`pwd`
	basename=`basename $file`
	echo $basename
	$2 -s $basename -no_optH false -ignore_unrecognized_res true -out:pdb >score_jd2.log 2>&1
	if [ $? -eq 0 ];then
		mv score.sc jd2_score.sc
		cd $pwd
		echo $fpwd/${basename%.pdb}_0001.pdb >> $1/jd2_pdb.list
	fi
	cd $pwd
done


