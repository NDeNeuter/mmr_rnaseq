import HTSeq
import os
import collections
import glob
import sys
import getopt
from multiprocessing import Pool
from subprocess import call

'''
HTSeq Counting

Input: 
  directory with samfiles mapped with HISAT2
  annotation file

Output:
  countfile readcounts.txt

'''


def readcount(myfile):
    almnt_file = HTSeq.SAM_Reader(myfile)
    counts = collections.Counter()
    for almnt in almnt_file:
        if not almnt.aligned:
            counts[ "_unmapped" ] += 1
            continue

        gene_ids = set()

        for iv, val in exons[ almnt.iv ].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[ gene_id ] += 1
        elif len(gene_ids) == 0:
            counts[ "_no_feature" ] += 1
        else:
            counts[ "_ambiguous" ] += 1

    return myfile,counts



## Command line Arguments
bamdirectory = ''
annotation_file = ''
outputfilename = ''

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"i:a:o:",["inputdir=","annotationfile=","outputfile="])
except getopt.GetoptError:
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--inputdir"):
            bamdirectory = arg
    elif opt in ("-o", "--outputfile"):
            outputfilename = arg
    elif opt in ("-a", "--annotationfile"):
            annotation_file = arg

#bamfile to be processed ends with:
filename_pattern= '.sam'

feature_type= ["gene","pseudogene"]
nthreads= 6

gtf_file = HTSeq.GFF_Reader(annotation_file)
exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)

for feature in gtf_file:
    if feature.type in feature_type:
        exons[feature.iv] += feature.attr["ID"]

os.chdir(bamdirectory)

listfiles = []

for mydir, subdirs, myfiles in os.walk(bamdirectory):
    for myfile in myfiles:
        if myfile.endswith('.sam'):
            #myfilepath = mydir.split('/')[-1]+'/'+myfile
            myfilepath = mydir + "/" + myfile
            listfiles.append(myfilepath)

print 'files to process:'

for myfile in listfiles:
    print '-' + myfile

print
print 'processing:'

p = Pool(nthreads)
summary = p.map(readcount, listfiles)
p.close()
p.join()

sumkeys=[]

print 'processing finished, writing to file...',

for j in summary:
    print j[0]
    for i in j[1]:
        if i not in sumkeys:
            sumkeys.append(i)



with open(outputfilename,'w') as countfile:
    countfile.write('genename\t')
    for file in summary:
        countfile.write(str(file[0])+ '\t')
    countfile.write('\n')

    for gene in sumkeys:
        countfile.write(str(gene) + '\t')
        for bamfile in summary:
            if gene in bamfile[1]:
                countfile.write(str(bamfile[1][gene])+'\t')
            else:
                countfile.write(str(0)+'\t')

        countfile.write('\n')

print 'finished'


