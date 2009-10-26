## The following configuration variables can be set in the file config.mak,
## otherwise the values in this file will be used as the defaults

include $(wildcard config.mak)

## ##############################
## DATA SET CONFIGURATION
##
COHORT_ID ?= 10
BIOINTDB ?= bioInt
PATHWAY ?= all
NULLITERS ?= 1000

## ##############################
## EXECUTABLE CONFIGURATION
##
DIGMA_INSTALL ?= /hive/users/${USER}/bin
KENT_PREFIX ?= ${HOME}

MAKECLUSTERFILES ?= ${KENT_PREFIX}/kent/src/hg/instinct/bioInt2/makeClusterFiles
MERGEEM ?= ${DIGMA_INSTALL}/helperScripts/mergeEMSwarmFiles.py
MAKEJOBSEM ?= ${DIGMA_INSTALL}/helperScripts/writeJobsEM.sh
MAKEJOBS ?= ${DIGMA_INSTALL}/helperScripts/writeJobs.sh
MERGE ?=${DIGMA_INSTALL}/helperScripts/mergeSwarmFiles.py

all: learningComplete inferenceComplete mergeFiles/done

mergeFiles/done: inferenceComplete
	mkdir -p mergeFiles
	/data/home/common/bin/python ${MERGE} outputFiles mergeFiles
	touch $@

inferenceComplete: ipara
	@echo 
	@echo "INSTRUCTIONS ---"
	@echo "INSTRUCTIONS --- Run the parasol queue in lpara via:"
	@echo "INSTRUCTIONS ---     para -batch=$< push"
	@echo "INSTRUCTIONS --- After the job is complete execute:"
	@echo "INSTRUCTIONS ---     touch $@"
	@echo "INSTRUCTIONS ---     make mergeFiles/done"

# parasol run for inference
ipara: jobs.list
	para -batch=$@ create $<

jobs.list: clusterFiles outputFiles/config.txt
	${JOBS} clusterFiles outputFiles > $@

outputFiles/config.txt: outputFiles outputFilesEM learningComplete
	${MERGEEM} outputFilesEM outputFiles

outputFiles:
	mkdir $@

learningComplete: lpara
	@echo 
	@echo 
	@echo 
	@echo "INSTRUCTIONS ---"
	@echo "INSTRUCTIONS --- Run the parasol queue in lpara via:"
	@echo "INSTRUCTIONS ---     para -batch=$< push"
	@echo "INSTRUCTIONS --- After the job is complete execute:"
	@echo "INSTRUCTIONS ---     touch $@"
	@echo "INSTRUCTIONS ---     make ipara"
	@false

lpara: jobsEM.list
	para -batch=$@ create $<

jobsEM.list: clusterFiles outputFilesEM
	${MAKEJOBSEM} $^ > $@

outputFilesEM:
	mkdir $@

clusterFiles:
	mkdir -p $@
	${MAKECLUSTERFILES} ${BIOINTDB} ${PATHWAY} $@ ${COHORT_ID} ${NULLITERS}
