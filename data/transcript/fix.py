


outFile = open('cdf_transcription.sub.format.csv','w')
with open('cdf_transcription.sub.tsv','r') as inFile:

	firstLine = inFile.readline().split()
	firstLine = ','.join(firstLine)
	outFile.write(firstLine)

	for line in inFile:
		line = line.split()
		gene = line[0].split('|')[0]
		entry = [gene]
		for x in line[1:]:
			entry.append(str(float(x)))
		entry = ','.join(entry) + '\n'
		outFile.write(entry)

outFile.close()