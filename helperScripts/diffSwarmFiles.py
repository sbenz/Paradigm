import sys, os, string, fnmatch

tolerance = 1e-5

def addData(inFile, sampleData):
    currentId = None
    for line in inFile:
        if line.startswith('>'):
            currentId = line[:-1].strip('>').strip()
            continue

        if currentId is None:
            continue

        data = line[:-1].split('\t')
        name = data[0]
        val = data[1]

        if not currentId in sampleData:
            sampleData[currentId] = {}

        sampleData[currentId][name] = val

    inFile.close()

def getInputStream(source):
    """Try to open a file, first using stdin if the source is the
    string -, second trying the file as a URL, third try as a regular
    file on a file system, and finally treat it just as a raw string.
    Based on toolbox.openAnything() at
    http://diveintopython.org/xml_processing/index.html"""
    if hasattr(source, "read"):
        return source
    import urllib
    try:
        return urllib.urlopen(source)
    except (IOError, OSError):
        pass
    try:                                  
        return open(source)               
    except (IOError, OSError):            
        pass                              
    import StringIO                       
    return StringIO.StringIO(str(source)) 

def main(fileA, fileB):
    dataA = {}
    addData(fileA, dataA)
    dataB = {}
    addData(fileB, dataB)
    
    if dataA.keys() != dataB.keys():
        print "Different sample names:"
        for s in set(dataA.keys()) ^ set(dataB.keys()): print s
    else:
        for s in dataA.keys():
            if (dataA[s].keys() != dataB[s].keys()):
                print "Entities differ for sample:", s
                for name in set(dataA[s].keys()) ^ set(dataB[s].keys()):
                    print ("\t%s" % name)
            else:
                outputSample = False
                for name, value in dataA[s].iteritems():
                    if (float(value) - float(dataB[s][name])) > tolerance:
                        if not outputSample:
                            print ">", s
                            outputSample = True
                        print "\t".join([s, value, dataB[s][name]])

def usage():
    print "python diffSwarmFiles.py file_a file_b"
    sys.exit(0)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        
    filea = getInputStream(sys.argv[1])
    fileb = getInputStream(sys.argv[2])

    main(filea, fileb)

    filea.close()
    fileb.close()
        
