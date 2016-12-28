# proCompare

Tool for comparing profile files, two by two, to analyse allele segregation. 

With the profile file and a group file, separating the isolates into groups, this scrip compares the loci on the profile with every two groups, indicating the number of loci that share alleles for the groups being compared. 

## Usage
	usage: proCompare.py [-h] [-p PROFILE] [-g GROUP]

	Profile comparison script.

	optional arguments:
		-h, --help				show this help message and exit
		-p PROFILE, --profile PROFILE
								Input profile tab file.
 		 -g GROUP, --group GROUP 
 		 						Group csv file
## Input

This scpit requires two input files:
- A tab separated file with the profile (MLST/cgMLST).
- A comma separated file, with isolates, identical in the two files, and the groups where they belong. The rows should correspond to the isolates and the columns to the different groups, using "0" to indicate absence and "1" to indicate presence of the group. The top left cell should be left blank.

It should look something like this:

|         | Group1 | Group2 | ... | GroupN |
| ------- | ------ | ------ | --- | ------ |
| Strain1 | 1      | 0      | ... | 1      |
| Strain2 | 1      | 1      | ... | 0      |
| Strain3 | 0      | 0      | ... | 1      |
| ...     | ...    | ...    | ... | ...    |
| StrainN | 1      | 0      | ... | 0      |

##Output

Comparing the groups two by two, it indicated the number of loci with shared alleles and the number of loci with only exclusive alleles for the two groups.


# proCompare_absent

Tool for comparing profile files, two by two, to analyse allele segregation. This script allows for missing data to be present in the profile. 

With the profile file and a group file, separating the isolates into groups, this scrip compares the loci in common on the profile with every two groups, indicating the number of loci that share alleles for the groups being compared. 

## Usage
	usage: proCompare.py [-h] [-p PROFILE] [-g GROUP]

	Profile comparison script.

	optional arguments:
		-h, --help				show this help message and exit
		-p PROFILE, --profile PROFILE
								Input profile tab file.
 		-g GROUP, --group GROUP 
 		 						Group csv file
 		-a ABSENT, --absent ABSENT
 								character representing th missing data

## Input

This scpit requires two input files:
- A tab separated file with the profile (MLST/cgMLST).
- A comma separated file, with isolates, identical in the two files, and the groups where they belong. The rows should correspond to the isolates and the columns to the different groups, using "0" to indicate absence and "1" to indicate presence of the group. The top left cell should be left blank.

It should look something like this:

|         | Group1 | Group2 | ... | GroupN |
| ------- | ------ | ------ | --- | ------ |
| Strain1 | 1      | 0      | ... | 1      |
| Strain2 | 1      | 1      | ... | 0      |
| Strain3 | 0      | 0      | ... | 1      |
| ...     | ...    | ...    | ... | ...    |
| StrainN | 1      | 0      | ... | 0      |

- The character that represents the missing data in the profile.

##Output

Comparing the groups two by two, it indicated the number of loci in common with shared alleles and the number of loci with only exclusive alleles for the two groups.

# Contact
Catarina Mendes (cimendes@medicina.ulisboa.pt)
