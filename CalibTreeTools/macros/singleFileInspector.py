from ROOT import *
import json, os, sys

print "Loading 'printSizes'"
gSystem.Load("printSizes_C.so")


def getEntry(fName):
  entry = {}
  entry['fName']   = fName
  entry['rNumber'] = fName.split("/")[-1].replace(".root","").split("_")[1]
  f = TFile(fName)
  entry['fSize']   = f.GetSize()

  treeList = []


  primaryKeys = f.GetListOfKeys()
  for key in primaryKeys:
    directory = f.Get(key.GetName())
    directoryNames = []  # Name of trees in directory

    for treeKey in directory.GetListOfKeys():
      if not treeKey.GetName() in directoryNames:
        directoryNames.append(treeKey.GetName())

        tree = directory.Get(treeKey.GetName())
        treeName = key.GetName()+"/"+treeKey.GetName()
        treeSize = sizeOnDisk(tree)

        treeEntry = {}
        treeEntry['Name']    = treeName
        treeEntry['Size']    = treeSize
        treeEntry['Entries'] = tree.GetEntries()
        branchList        = []

        for branchKey in tree.GetListOfBranches():
          branchName = branchKey.GetName()
          branch     = tree.GetBranch(branchName)
          branchSize = sizeOnDisk(branch,1)

          branchEntry = {}
          branchEntry["Name"] = branchName
          branchEntry["Size"] = branchSize
          branchEntry["Entries"] = branch.GetEntries()
          branchList.append(branchEntry)

        treeEntry["Branches"] = branchList
        treeList.append(treeEntry)
  entry["Trees"] = treeList
  return entry


if len(sys.argv) < 3:
  print "Error, expecting 2 arguments (input and output files)!"
  sys.exit()
outputFile     = sys.argv[2]
inputFile      = sys.argv[1]

print "%s %s"%(outputFile,inputFile)

jsonData = [getEntry(inputFile)]
print jsonData[0]


with open(outputFile,'w') as data_file:
  data_file.write(json.dumps(jsonData))
