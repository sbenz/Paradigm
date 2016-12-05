#!/usr/bin/python
import sys, os, string, fnmatch, re

skipEvidenceOptions = ["epsilon", "epsilon0", "reverse"]

def readEvidenceConfiguration(filename):
    f = open(filename, "r")
    typeToOpts = {}
    infLines = []
    for x in f.readlines():
        m = re.match("(\w+)\s+\[(.*)\]", x)
        if m == None: sys.exit("Bad configuration line:\n %s" % x)
        if m.group(1) == "inference": infLines.append(x)
        elif m.group(1) == "evidence":
            opts = {}
            for f in m.group(2).split(","):
                ff = f.split("=", 1)
                opts[ff[0]] = ff[1]
            typeToOpts[opts["suffix"]] = opts
    return infLines, typeToOpts

def getFilesMatching(baseDir, patterns):
    list = []
    
    for root, dirs, files in os.walk(baseDir):
        for file in files:
            ptr = os.path.join(root, file)        
            for pattern in patterns:
                if fnmatch.fnmatch(ptr, pattern):
                    list.append(ptr)
                    
    return list 


def addData(fname, sampleData):
    inFile = open(fname)

    inFile.readline()
    type = None
    for line in inFile:
        if line.startswith('>'):
            type = line.split("'").pop(1)
            continue

        if type is None:
            continue

        data = line.strip("\n").split('\t')
        parent = int(data[1])
        child = int(data[0])
        val = float(data[2])

        if type not in sampleData:
            sampleData[type] = []
            for i in range(3):
                sampleData[type].append([0,0,0])

        sampleData[type][parent][child] += val

    inFile.close()
    
def mergeFiles(outdirectory, files, infLines, evidLines):

    sampleData = {}
    for f in files:
        addData(f, sampleData)

    outname = os.path.join(outdirectory, "config.txt")
    outfile = open(outname,'w')
    for l in infLines: outfile.write(l)

    for type in sampleData:
        evid = evidLines[type]
        outfile.write("evidence [")
        for k,v in evid.items():
            if k not in skipEvidenceOptions: outfile.write("%s=%s,"% (k,v))
        outfile.write("factorParams=")
        for i in range(3):
            for j in range(3):
                sampleData[type][i][j] /= len(files)
                if i>0 or j>0:
                    outfile.write(";")
                outfile.write(str(sampleData[type][i][j]))
        outfile.write("]\n")
    
def main(indirectory, outdirectory):
    infLines, evidConf = readEvidenceConfiguration("%s/config.txt" % indirectory)
    allfiles = getFilesMatching(indirectory, ["*learned_parameters.fa"])
    print "found ", len(allfiles), " files total"
    
    mergeFiles(outdirectory, allfiles, infLines, evidConf)
    

def usage():
    print "python mergeSwarmFiles.py indirectory outdirectory"
    sys.exit(0)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        
    indirectory = sys.argv[1]
    outdirectory = sys.argv[2]

    main(indirectory, outdirectory)
        
