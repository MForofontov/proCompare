# proCompare

Tool for comparing profile files, two by two, to analyse allele segregation. 

With the profile file and a group file, separating the isolates into groups, this scrip compares the loci on the profile with every two groups, indicating the number of loci that share alleles for the groups being compared. 

## Usage
	usage: python proCompare.py -h 
usage: proCompare.py [-h] -p PROFILE -g GROUP -o OUTDIR -a ABSENT

Created on Wed Apr 20 17:05:00 2022

@author: Inês Mendes, Mykyta Forofontov

optional arguments:
  -h, --help            show this help message and exit
  -p PROFILE, --profile PROFILE
                        Input profile tab file
  -g GROUP, --group GROUP
                        group tsv file
  -o OUTDIR, --outdir OUTDIR
                        outdir path
  -a ABSENT, --absent ABSENT
                        character representing missing data
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

## Output

Comparing the groups two by two, it indicated the number of loci with shared alleles and the number of loci with only exclusive alleles for the two groups.

# Contact
Inês Mendes (cimendes@medicina.ulisboa.pt)
