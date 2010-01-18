#!/usr/bin/env bash
#for BASH shells

if [ $# -lt 2 ] ; then
    echo "Usage: "
    echo "  writeJobs.sh indirectory outdirectory"
    exit 1
fi

TMP=$1
DIRECTORY=`echo $TMP | sed "s/\///g"`

OUTDIR=$2

# The loop will get all files with .fa extension in the /files directory
for P in ${DIRECTORY}/*_pathway.tab; do
    INFIX=`echo $P | sed "s/_pathway.tab//g" | sed "s/${DIRECTORY}\///g"`
    
    for F in ${DIRECTORY}/*${INFIX}_*CNV.tab; do
	T=`echo $F | sed "s/_CNV.tab//g"`
	B=`basename $T`
	echo "/hive/users/${USER}/bin/paradigm -p ${P} -b ${T} -c ${OUTDIR}/config.txt -o $OUTDIR/${B}_output.fa"
    done

done
