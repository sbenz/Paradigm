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


def groupFilesById(allfiles):

    groups = {}
    for file in allfiles:
        i = file.find("pid")
        if i == -1:
            continue

        s = file[i:].split('_')
        pid = s[0] + '_' + s[1]

        if not pid in groups:
            groups[pid] = []

        groups[pid].append(file)

    return groups

def addData(fname, sampleData):
    inFile = open(fname)

    prefix = ""
    if fname.find("nw_") != -1:
        prefix = "nw_"
    elif fname.find("na_") != -1:
        prefix = "na_"
        
    currentId = None
    for line in inFile:
        if line.startswith('>'):
            currentId = prefix + line[:-1].strip('>').strip()
            continue

        if currentId is None:
            continue

        data = line[:-1].split('\t')
        name = data[0]
        val = data[1]

        if not name in sampleData:
            sampleData[name] = {}

        sampleData[name][currentId] = val

    inFile.close()
    
def outputData(outname, sampleData):
    outFile = open(outname, 'w')

    entities = sampleData.keys()
    entities.sort()
    
    s = "id" + '\t' + '\t'.join(entities) + '\n'
    outFile.write(s)
    
    sampleNames = []
    for e in entities:
        names = sampleData[e].keys()
        for n in names:
            if not n in sampleNames:
                sampleNames.append(n)

    sampleNames.sort()
    sampleNames.reverse() # put "sample" on top

    for sample in sampleNames:
        data = []
        for e in entities:
            try:
                val = sampleData[e][sample]
            except:
                val = "NA"
            data.append(val)
        s = sample + '\t' + '\t'.join(data) + '\n'
        outFile.write(s)

    outFile.close()

def outputDataTranspose(outname, sampleData):
    outFile = open(outname, 'w')

    entities = sampleData.keys()
    entities.sort()
    
    sampleNames = []
    for e in entities:
        names = sampleData[e].keys()
        for n in names:
            if not n in sampleNames:
                sampleNames.append(n)

    sampleNames.sort()
    sampleNames.reverse() # put "sample" on top

    s = "id" + '\t' + '\t'.join(sampleNames) + '\n'
    outFile.write(s)
    
    for e in entities:
        data = []
        for sample in sampleNames:
            try:
                val = sampleData[e][sample]
            except:
                val = "NA"
            data.append(val)
        s = e + '\t' + '\t'.join(data) + '\n'
        outFile.write(s)

    outFile.close()
    
    
def mergeGroups(outdirectory, groups):

    for g,files in groups.iteritems():
        sampleData = {}
        for f in files:
            addData(f, sampleData)

        outname = os.path.join(outdirectory, "merged_" + g + ".out")

        print "merging", len(files), "files into", outname
        outputData(outname, sampleData)

        outnameTranspose = os.path.join(outdirectory, "merged_transpose_" + g + ".out")
        outputDataTranspose(outnameTranspose, sampleData)

def main(indirectory, outdirectory):
    allfiles = getFilesMatching(indirectory, ["*.fa"])
    print "found ", len(allfiles), " files total"
    
    groups = groupFilesById(allfiles)
    print "grouped files into ", len(groups), " groups"

    mergeGroups(outdirectory, groups)
    

def usage():
    print "python mergeSwarmFiles.py indirectory outdirectory"
    sys.exit(0)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        
    indirectory = sys.argv[1]
    outdirectory = sys.argv[2]

    main(indirectory, outdirectory)
        
