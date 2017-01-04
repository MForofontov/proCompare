import argparse, csv, sys, numpy, time, os
from Bio import SeqIO
import itertools

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

	group_dic={}
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

			group_dic[name_trait] = p.keys()
	
	print "Comparing groups:"
	for key, value in group_dic.items():
		print key
		print " ...containing %s isolates" % str(len(value))
	print ''

	groups=list(group_dic.keys())

	return group_dic, groups

def loadProfile(filename,dataStructure):
	'''Load profile file into the program, separating the samples into the groups'''

	with open(filename, 'r') as profile:
		reader = csv.reader(profile, delimiter='\t')
		genes=reader.next()[1:]
		profileDic={}

		for line in reader:
			skipIsolate=True
			strain=line[0]
			profile=line[1:] #has the same order as in genes

			for key, value in dataStructure.items():
				if strain in value:
					group=key
					skipIsolate=False
					break
				else:
					pass

			if not skipIsolate:
				for i in range (0,len(genes)):
					if genes[i] not in profileDic.keys():
						temp={}
						temp[group]=[profile[i]]
						profileDic[genes[i]]=temp
					else:
						temp=profileDic[genes[i]]
						if group in temp.keys():
							temp[group].append(profile[i])
						else:
							temp[group]=[profile[i]]
					
						profileDic[genes[i]]=temp
			else:
				pass

	print "profile size: " + str(len(profileDic))

	return profileDic
	
def proCompare(dataStructure, key1, key2, absent):
	'''Comparison of the the two groups in the profile'''

	shared_allele=[]
	exclusive_alleles=[]

	#performing counts
	for gene, value in dataStructure.items():

		if absent in dataStructure[gene][key1] or absent in dataStructure[gene][key2]:
			pass
		else:

			#group one
			group_one=set(value[key1])
			#group two
			group_two=set(value[key2])

			#common alleles in the two groups
			common_alleles=len(group_one.intersection(group_two))
			if common_alleles != 0:
				shared_allele.append(gene)
			else:
				exclusive_alleles.append(gene)

			#Alleles present only in group one
			group_one_diff=len(group_one.difference(group_two))

			#Alleles present only in group two
			group_two_diff=len(group_two.difference(group_one))

			value['Stats']=[group_one,group_two,common_alleles,group_one_diff,group_two_diff]

	print '\n Comparing: ' + str(key1) + ' + ' + str(key2) + ' with ' + str(len(shared_allele)+len(exclusive_alleles)) + ' loci in common'
	print "\tLoci with shared alleles: " + str(len(shared_allele))
	print "\tLoci with exclusive alleles: " + str(len(exclusive_alleles))

	return dataStructure

def main():

	parser = argparse.ArgumentParser(description='Profile comparison script, with missing data.', epilog='by Catarina Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-p', '--profile', help='Input profile tab file.')
	parser.add_argument('-g', '--group', help='group csv file')
	parser.add_argument('-a', '--absent', help='character representing missing data')

	args = parser.parse_args()

	sys.stdout = Logger("./")
	#Loads group file
	groups_dic, groups_list=loadGroups(args.group)

	#Loads profile file
	profiles=loadProfile(args.profile,groups_dic)

	#performs profile comparison for two groups at a time
	for a, b in itertools.combinations(groups_list, 2):
		genes=proCompare(profiles, a, b, args.absent)

	print "\nFinished"
	sys.exit(0)
	
if __name__ == "__main__":
    main()