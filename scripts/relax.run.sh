#!/bin/bash
#inputs:
#	$1	pdb dictionary
#	$2  mpi cores (if you want to disable mpi ,input arg should be 0)
#   $3  flag file
#   $4  exec file
#   $5  jd_score.sc.list (which contain the path of score_jd's score.sc)
#   $6  pdb dictionary
#outputs:
#	{PDB DIR}/{PDB}_{index}.pdb
#	{PDB DIR}/relax_sc.list

pwd=`pwd`
tmpfile=$$.fifo
mkfifo $tmpfile
exec 4<>$tmpfile
rm $tmpfile

mpi_cores=$1
threads=$2
relax_flag=$3
exec_file=$4
jd2_file=$5
pdb_dir=$6
# relax_flag=$(dirname $(readlink -f "$0"))/relax.flag
{
for (( i=1;i<=${threads};i++ ))
do
	echo;
done
} >&4
#for file in `find $1 -name "*.pdb" |grep -v "init"`
if [ -d jd2_file ];then
find $jd2_file -name "*.pdb" > $jd2_file/jd2_pdb.list
file_path=$jd2_file/jd2_pdb.list
else
file_path=$jd2_file
fi
for file in `cat $file_path`
do
	read
	dir=`dirname $file`
	cd $dir
	basename=`basename $file`
	echo $basename
	if [ $mpi_cores -eq 0 ]
	then
		($exec_file -s $basename @ $relax_flag >relax.log 2>&1 </dev/null;echo >&4) &
	else
		(mpirun -np $mpi_cores --use-hwthread-cpus $exec_file -s $basename @ $relax_flag < /dev/null;echo >&4) &
	fi
	cd $pwd

done <&4
wait 
exec 4>&-
# find $1 -name "score.sc" |xargs -i mv {} relax_scores.sc
find $pdb_dir -name "score.sc" > $pdb_dir/relax_sc.list

