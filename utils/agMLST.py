import argparse, csv, sys

def main():

	parser = argparse.ArgumentParser(description='agMLST profile file from cg and wg MLST profile files.', epilog='by Catarina Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-wg', '--wgMLST', help='wgMLST profile tab file.')
	parser.add_argument('-cg', '--cgMLST', help='cgMLST profile tab file.')

	args = parser.parse_args()

	core_genes=[]
	with open(args.cgMLST, 'r') as cgMLST:
		reader = csv.reader(cgMLST, delimiter='\t')
		core_genes=reader.next()

	all_genes=[]
	with open(args.wgMLST, 'r') as wgMLST:
		reader=csv.reader(wgMLST, delimiter='\t')
		all_genes=reader.next()

	with open(args.wgMLST, 'r') as wgMLST:
		with open('agMLST.tsv', 'w') as agMLST:
			reader=csv.reader(wgMLST, delimiter='\t')
			for line in reader:
				toWrite=[line[0]]
				for i in range(0,len(line)):
					gene=all_genes[i]
					if gene in core_genes:
						pass
					else:
						toWrite.append(line[i])

				agMLST.write('\t'.join(toWrite) + '\n')

	print "\nFinished"
	sys.exit(0)
	
if __name__ == "__main__":
    main()