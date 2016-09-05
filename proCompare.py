import argparse, csv, sys, numpy, time, os
from Bio import SeqIO

class Logger(object):
	'''Class to create a log file of all that is printed onto the console'''

	def __init__(self, out_directory):
		self.logfile = os.path.join(out_directory, "run.log")
		if os.path.isfile(self.logfile):
			print "Logfile already exists! It will be overwritten..." + "\n"
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")
	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()
	def flush(self):
		pass


def loadGroups(filename):
	'''Load group file into a dictionary'''

	r={}
	with open(filename, 'r') as groups:
		reader = zip(*csv.reader(groups, delimiter=','))

		for group in range(1,len(reader)):
			p=dict(zip(reader[0], reader[group]))
			if "" in p:
				name_trait = p[""]
				del p[""]
			elif "Name" in p:
				name_trait = p["name"]
				del p["Name"]
			else:
				sys.exit("Make sure the top-left cell in the traits file is either empty or 'Name'.")

			for key,value in p.items():
				if value == '0':
					del p[key]

			r[name_trait] = p.keys()
	
	print "Comparing groups:"
	for key, value in r.items():
		print key
		print " ...containing %s isolates" % str(len(value))
	print ''

	return r

def loadProfile(filename,dataStructure):
	'''Load profile file into the program, separating the samples into the groups'''

	with open(filename, 'r') as profile:
		reader = csv.reader(profile, delimiter='\t')
		genes=reader.next()[1:]
		geneDic={}

		for line in reader:
			skipIsolate=False
			strain=line[0]
			profile=line[1:] #has the same order as in genes

			for key, value in dataStructure.items():
				if strain in value:
					group=key
				else:
					skipIsolate=True

			if skipIsolate:
				for i in range (1,len(genes)):
					if genes[i] not in geneDic.keys():
						temp={}
						temp[group]=[profile[i]]
						geneDic[genes[i]]=temp
					else:
						temp=geneDic[genes[i]]
						if group in temp.keys():
							temp[group].append(profile[i])
						else:
							temp[group]=[profile[i]]
					
						geneDic[genes[i]]=temp
			else:
				pass


	return geneDic

def printCSV(genes,numAleles,alleles_G1,alleles_G2,common,exclusive_G1,exclusive_G2,diff_Alleles):
	'''Method to print a TSV with the profile comparisons for the two groups. This file can be easily manipulated in MS Excel'''
	
	fh=open("results.tsv", 'w')

	fh.write('Stats'+'\t'+"\t".join(genes)+'\n')
	fh.write('Unique Alleles'+'\t'+'\t'.join(numAleles)+'\n')
	fh.write('Unique Horse Alleles'+'\t'+'\t'.join(alleles_G1)+'\n')
	fh.write('Unique Human Alleles'+'\t'+'\t'.join(alleles_G2)+'\n')
	fh.write('Unique Common Alleles'+'\t'+'\t'.join(common)+'\n')
	fh.write('Unique Exclusive Horse Alleles'+'\t'+'\t'.join(exclusive_G1)+'\n')
	fh.write('Unique Exclusive Human Alleles'+'\t'+'\t'.join(exclusive_G2)+'\n')
	str_diffAlleles=[]
	for item in diff_Alleles:
		str_diffAlleles.append(str(item))
	fh.write('Sum Diff Alleles'+'\t'+'\t'.join(str_diffAlleles))

	fh.close()
	
def proCompare(dataStructure):
	'''Comparison of the the two groups in the profile'''

	genes=[]
	numAleles=[]
	alleles_G1=[]
	alleles_G2=[]
	common=[]
	exclusive_G1=[]
	exclusive_G2=[]
	diff_Alleles=[]

	shared_allele=[]
	exclusive_alleles=[]

	#performing counts
	for gene, value in dataStructure.items():

		genes.append(gene)

		groups = value.keys()

		#group one
		group_one=set(value[groups[0]])
		alleles_G1.append(str(len(group_one)))

		#group two
		group_two=set(value[groups[1]])
		alleles_G2.append(str(len(group_two)))


		allAlleles=len(group_one.union(group_two))
		numAleles.append(str(allAlleles))

		#common alleles in the two groups
		common_alleles=len(group_one.intersection(group_two))
		common.append(str(common_alleles))
		if common_alleles != 0:
			shared_allele.append(gene)
		else:
			exclusive_alleles.append(gene)

		#Alleles present only in group one
		group_one_diff=len(group_one.difference(group_two))
		exclusive_G1.append(str(group_one_diff))

		#Alleles present only in group two
		group_two_diff=len(group_two.difference(group_one))
		exclusive_G2.append(str(group_two_diff))

		#Sum of exclusive alleles for each group
		diff_Alleles.append(int(group_one_diff + group_two_diff))

		value['Stats']=[group_one,group_two,common_alleles,group_one_diff,group_two_diff]

	printCSV(genes,numAleles,alleles_G1,alleles_G2,common,exclusive_G1,exclusive_G2,diff_Alleles)

	print "profile size: " + str(len(dataStructure))
	print ''
	print "loci with shared alleles within the two groups: " + str(len(shared_allele))
	print "loci with exclusive alleles for the two groups: " + str(len(exclusive_alleles))
	print ''

	return dataStructure


def checkGenome(dataStructure, schemaDir):
	'''Checks the size of the genome with loci that have shared alleles in the two groups'''

	schemaFiles=os.listdir(schemaDir)

	gcGenes=[]
	cgSizes=[]
	commonSizes=[]

	for gene in schemaFiles:
		if gene in dataStructure.keys():
			gcGenes.append(gene)

	for f in gcGenes:

		with open(schemaDir+'/'+f, "r") as file:
			seq=SeqIO.read(file,"fasta")

			#all core
			cgSizes.append(int(len(seq)))

			#common loci
			stats = dataStructure[f]['Stats']
			if stats[2] != 0:
				commonSizes.append(int(len(seq)))

	print "Total size of core genome: " + str(sum(cgSizes))+' bp'
	print "Medium size of loci in the core genome: " + str(numpy.median(numpy.array(cgSizes)))+' bp'
	print ''
	print "Total size of loci with common alleles in the two groups:: " + str(sum(commonSizes))+' bp'
	print "Medium size of loci with common alleles in the two groups: " + str(numpy.median(numpy.array(commonSizes)))+' bp'


def main():

	parser = argparse.ArgumentParser(description='Profile comparison script.', epilog='by Catarina Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-p', '--profile', help='Input profile tab file.')
	parser.add_argument('-g', '--group', help='group csv file')
	parser.add_argument('-s', '--schema', help='directory for the schema')


	args = parser.parse_args()

	sys.stdout = Logger("./")

	#Loads group file
	groups=loadGroups(args.group)

	#Loads profile file
	profiles=loadProfile(args.profile,groups)

	#performs profile comparison for the two groups in the group file. Prints TSV report.
	genes=proCompare(profiles)

	#checks total size of the loci with shared alleles within the two groups
	checkGenome(genes,args.schema)

	print "\nFinished"
	sys.exit(0)
	
if __name__ == "__main__":
    main()