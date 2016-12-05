#!/data/home/common/bin/python

import os
import sys
import fnmatch
import sqlalchemy

resultMatrix = {} # first key is entity, then sample
sampleList = []

def usage():
	print "Usage: "+sys.argv[0]+" merged_dir"
	sys.exit(0)

def getFilesMatching(baseDir, patterns):
	list = []
	
	for root, dirs, files in os.walk(baseDir):
		for file in files:
			ptr = os.path.join(root, file)
			for pattern in patterns:
				if fnmatch.fnmatch(ptr, pattern):
					list.append(ptr)
	return list

def addFileToResult(file):
	pid = file[:-4].split("_").pop()
	fh = open(file,"r")
	header = fh.readline().strip("\n").split("\t")
	sampleOrder = header[1:]
	for sample in sampleOrder:
		if sample.startswith("sample") and sample.split(" ").pop(0) not in sampleList:
			sampleList.append(sample.split(" ").pop(0))
	for line in fh:
		dataA = line.strip("\n").split("\t")
		entity = dataA.pop(0)
		resultMatrix[pid+"_"+entity] = {}
		for i in range(len(dataA)):
			if sampleOrder[i].startswith("sample"):
				resultMatrix[pid+"_"+entity][sampleOrder[i].split(" ").pop(0)] = dataA[i]
	fh.close()

def main(directory,outputFile):
	files = getFilesMatching(directory, ["*_pid_*"])

	print "Loading data",
	sys.stdout.flush()
	for f in files:
		pid = f[:-4].split("_").pop()
		if pid == "example":
			continue
		addFileToResult(f)
		print ".",
		sys.stdout.flush()
		
	print "Done reading data, converting sample ids...",
	sys.stdout.flush()	
	# open the sql connection
	engine = sqlalchemy.create_engine('mysql://tcgacat:tcga$TCGA@localhost/bioIntTCGAOV')
	connection = engine.connect()
	transSampleList = []
	for i in range(len(sampleList)):
		sampleID = sampleList[i].split("_").pop()
		query = "select name from samples where id = "+sampleID
		result = connection.execute(query)
		row = result.fetchone()
		transSampleList.append(row['name'])
	# close the sql connections cause we're done now
	connection.close()
	print "done."
	print "Printing Results..."
	resultFile = open(outputFile,"w")
	resultFile.write("pid_entity\t")
	resultFile.write("\t".join(transSampleList)+"\n")
	for entity in resultMatrix:
		resultFile.write(entity)
		for sample in sampleList:
			resultFile.write("\t")
			#print sample
			#print resultMatrix[entity]
			#sys.exit(0)
			if sample in resultMatrix[entity]:
				resultFile.write(resultMatrix[entity][sample])
		resultFile.write("\n")
				
	resultFile.close()
	

if __name__ == "__main__":
	if len(sys.argv) != 3:
		usage()

	directory = sys.argv[1]
	outputFile = sys.argv[2]
	main(directory,outputFile)

