import argparse, csv, sys, numpy, time, os
from Bio import SeqIO

class Logger(object):
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
	for key in r.keys():
		print key
	print ''

	return r

def loadProfile(filename,dataStructure):

	with open(filename, 'r') as profile:
		reader = csv.reader(profile, delimiter='\t')
		genes=reader.next()[1:]
		geneDic={}

		for line in reader:
			strain=line[0]
			profile=line[1:] #has the same order as in genes

			for key, value in dataStructure.items():
				if strain in value:
					group=key

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

	return geneDic

def printCSV(genes,numAleles,alleles_G1,alleles_G2,common,exclusive_G1,exclusive_G2,diff_Alleles):
	
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
		#print gene
		groups = value.keys()

		#allAlleles=set(set(value[groups[0]]),set(value[groups[1]]))

		group_one=set(value[groups[0]])
		alleles_G1.append(str(len(group_one)))
		#print 'Alleles ' +  str(groups[0]) + ': ' + str(len(group_one))
		group_two=set(value[groups[1]])
		alleles_G2.append(str(len(group_two)))
		#print 'Alleles ' +  str(groups[1]) + ': ' + str(len(group_two))

		allAlleles=len(group_one.union(group_two))
		numAleles.append(str(allAlleles))
		#print 'All alleles: ' + str(allAlleles)

		common_alleles=len(group_one.intersection(group_two))
		common.append(str(common_alleles))
		#print 'Common Alleles: ' + str(common_alleles)
		if common_alleles != 0:
			shared_allele.append(gene)
		else:
			exclusive_alleles.append(gene)

		group_one_diff=len(group_one.difference(group_two))
		exclusive_G1.append(str(group_one_diff))
		#print 'Eclusive ' +  str(groups[0]) + ': '+ str(group_one_diff)

		group_two_diff=len(group_two.difference(group_one))
		exclusive_G2.append(str(group_two_diff))
		#print 'Eclusive ' +  str(groups[1]) + ': '+ str(group_two_diff)

		value['Stats']=[group_one,group_two,common_alleles,group_one_diff,group_two_diff]

		diff_Alleles.append(int(group_one_diff + group_two_diff))
	
	printCSV(genes,numAleles,alleles_G1,alleles_G2,common,exclusive_G1,exclusive_G2,diff_Alleles)

	print "profile size: " + str(len(dataStructure))
	print ''
	print "loci with shared alleles within the two groups: " + str(len(shared_allele))
	print "loci with exclusive alleles for the two groups: " + str(len(exclusive_alleles))
	print ''
	#print "maximum number of different alleles per locus for the two groups: " + str(max(diff_Alleles))
	#print "minimum number of different alleles per locus for the two groups: " + str(min(diff_Alleles))
	#print ''
	#print "average number of different alleles per locus for the two groups:"
	#lala=sum(diff_Alleles) / float(len(diff_Alleles))
	#print lala

	return dataStructure


def checkGenome(dataStructure, schemaDir):

	schemaFiles=os.listdir(schemaDir)
	#print schemaFiles
	gcGenes=[]
	cgSizes=[]
	commonSizes=[]
	#count = 0
	for gene in schemaFiles:
		if gene in dataStructure.keys():
			gcGenes.append(gene)
	#print len(gcGenes) #994

	for f in gcGenes:

		with open(schemaDir+'/'+f, "r") as file:
			seq=SeqIO.parse(file,"fasta")

			for item in seq:
				#all core
				cgSizes.append(int(len(item)))

				#common loci
				stats = dataStructure[f]['Stats']
				if stats[2] != 0:
					commonSizes.append(int(len(item)))

	#print "Number of loci in core genome: " + str(len(cgSizes))
	#print "Number of loci with common alleles in the two groups: " + str(len(commonSizes))
	#print ''
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
	#parser.add_argument('--version', help='Display version, and exit.', default=False,action='store_true')

	args = parser.parse_args()

	start_time = time.time()

	sys.stdout = Logger("./")

	groups=loadGroups(args.group)
	#print groups

	profiles=loadProfile(args.profile,groups)
	#print profiles
	#print geneNames

	genes=proCompare(profiles)

	checkGenome(genes,args.schema)


	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken,3600)
	minutes, seconds = divmod(rest, 60)
	#print "Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

	print "\nFinished"
	sys.exit(0)
	
if __name__ == "__main__":
    main()