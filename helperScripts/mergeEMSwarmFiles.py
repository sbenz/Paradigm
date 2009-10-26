import sys, os, string, fnmatch

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
    
def mergeFiles(outdirectory, files):

    sampleData = {}
    for f in files:
        addData(f, sampleData)

    outname = os.path.join(outdirectory, "config.txt")
    outfile = open(outname,'w')
    outfile.write("inference [method=JTREE,updates=HUGIN,verbose=1]\n")

    for type in sampleData:
        node = type.split(".").pop(0)[1:]
        outfile.write("evidence [suffix="+type+",node="+node+",disc=-1.3;1.3,factorParams=")
        for i in range(3):
            for j in range(3):
                sampleData[type][i][j] /= len(files)
                if i>0 or j>0:
                    outfile.write(";")
                outfile.write(str(sampleData[type][i][j]))
        outfile.write("]\n")
    
def main(indirectory, outdirectory):
    allfiles = getFilesMatching(indirectory, ["*learned_parameters.fa"])
    print "found ", len(allfiles), " files total"
    
    mergeFiles(outdirectory, allfiles)
    

def usage():
    print "python mergeSwarmFiles.py indirectory outdirectory"
    sys.exit(0)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        
    indirectory = sys.argv[1]
    outdirectory = sys.argv[2]

    main(indirectory, outdirectory)
        
