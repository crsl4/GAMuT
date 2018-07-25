#!/bin/sh
# Tell the SGE that this is an array job, with "tasks" to be numbered 1 to 10000
#$ -t 1-1000
# When a single command in the array job is sent to a compute node,
# its task number is stored in the variable SGE_TASK_ID,
# so we can use the value of that variable to get the results we want.
# It needs to have run simulation-pre.r ahead of time

module load julia
module load R

PRJDIR="${HOME}/gamut"
#DATADIR="${PRJDIR}/data"
#OUTDIR="${PRJDIR}/output/FastQC"

TMP_NAME=`/usr/bin/mktemp -u XXXXXXXX`
TMPDIR="/scratch/${JOB_ID}_${SGE_TASK_ID}_${TMP_NAME}"
if [ -e $TMPDIR ]; then
   echo "Error. Scratch dir already exists on `hostname`"
   exit
else
    mkdir $TMPDIR
fi


cp ${PRJDIR}/*.jl ${TMPDIR}
cp ${PRJDIR}/*.r ${TMPDIR}
cd ${TMPDIR}

julia ${TMPDIR}/simulation-one.jl $SGE_TASK_ID

/bin/rm -f ${TMPDIR}/*.r
/bin/rm -f ${TMPDIR}/*.sh
/bin/rm -f ${TMPDIR}/*.jl


rsync -av ${TMPDIR}/ ${PRJDIR}

/bin/rm -fr ${TMPDIR}

module unload R
module unload julia 
