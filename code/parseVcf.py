# this code parses the ExAC aggregation vcf data file and extracts only the following columns:
# chr numer, position, AC, AC_AFR, AC_AMR, AC_Adj, AC_EAS, AC_FIN, AC_Het, AC_Hom, AC_NFE, AC_SAS, AF, 
# AN, AN_AFR, AN_AMR, AN_Adj, AN_EAS, AN_FIN, AN_NFE, AN_SAS
# for positions that pass the filtering
# the code outputs 22 samller files named chr.#_analysis

#!/usr/bin/env python

import sys
import os

inputFile = '~/release0.1/ExAC.r0.1.sites.vep.vcf'
outputDir = '~/data/filtered/'
chromosomesToAnalyze = range(1,23)
outSkeleton = "chr.%s_analysis"

def parseFile (chrNumber):

	outputFile = os.path.join (outputDir, outSkeleton % chrNumber)
	out = open(outputFile, 'w')
	ifs = open (inputFile)

	# go through
	for line in ifs.xreadlines():

		# drop everything that has a '#'
		if (line.startswith('#')):
			continue

		cols = line.split()

		if (cols[0] == 'X'):
			continue
			
		# check the index (chromosome)
		thisIdx = int(cols[0])

		if (thisIdx < chrNumber):
			continue

		if (thisIdx == chrNumber):
			# check if passes filtering
			if (cols[6] != 'PASS'):
				continue
			else:
				info = cols[7].split(';')
				out.write(cols[0] + ' ' + cols[1] + ' ' + ' '.join(info[0:19]))
				out.write('\n')

		if (thisIdx > chrNumber):
			break		

	ifs.close()		

def main():
	# go through chromosomes
	for c in chromosomesToAnalyze:
		parseFile(c)



if __name__ == "__main__":
	main()




