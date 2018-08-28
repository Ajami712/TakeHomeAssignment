import sys

NumArguments = list(sys.argv)

def Error_NoArguments():
	print '''Please supply the following arguments: 
1. Input_VCF_Directory
2. Output_Directory
3. (Optional) Set Verbose to "True" or "False" (boolean, default False)
4. (Optional) Set number of threads for API retrieval (integer, default 200)

Or provide the argument "help" for further reading.
	'''
def Help_Text():
	print '''#=======================================
VCF_Annotation
#=======================================

Input: A directory to a vcf file (Input_VCF)
Output: An annotated vcf written to a specified directory (Output_Filename)

The vcf is annotated with the following information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.)  - represented by TYPE in vcf file's INFO column
2. Depth of sequence coverage at the site of variation - represented by DP category in vcf file's INFO column
3. Number of reads supporting the variant - represented by AO category in vcf file's INFO column
4. Percentage of reads supporting the variant versus those supporting reference reads - calculated with value 3 / value 2
5. Allele frequency of variant from Broad Institute ExAC Project API (http://exac.hms.harvard.edu/)

Annotations are added to VCF files as individual columns.

Users have the option to set the parameter Verbose to True or False (default False).
When Verbose = True, all information in the original vcf file is retained in the annotation.
When Verbose = False, only the chromosome, position, reference allele, alternative allele, and quality values are retained.

The VCF_Annotation class object stores four objects: Metadata, ColNames, TextBody, and Indices.
Metadata contains the metadata of the input VCF file.
ColNames contains the header of the input VCF file.
TextBody contains the rows with variant information from the VCF file.
Indices defines the variant IDs used for the ExAC Project API, and is used to 
define associations between TextBody and information mined from the ExAC API.

Variants that have more than one alternative allele are decomposed to multiple lines.
According to http://samtools.github.io/hts-specs/VCFv4.1.pdf, "it is permitted to have multiple records with the same POS."

#=======================================
Important Functions within VCF_Annotation:
#=======================================

Run_Full_Annotation(Input_VCF, Output_Filename, Verbose = False, numThreads = 200):
Runs the full annotation pipeline by running Parse_VCF and Finalize_VCF.
The quick one-line-of-code way to run this analysis.
numThreads is the number of threads for communicating with the ExAC Project API - see Fast_ExAC_API_Mining.

Parse_VCF(Input_VCF):
Populates the Indices, TextBody, and ColNames data structures of the class object.

Finalize_VCF(Output_Filename, Verbose = False, numThreads = 200):
Runs Fast_ExAC_API_Mining, and integrates API results with VCF information.
This is done by converting TextBody to a pandas dataframe with colnames defined by ColNames and indices defined by Indices.
The API results are returned to the function as a hash table with Indices as keys.
If verbosity is False, unnecessary columns are dropped from the dataframe.
The final annotated file is written to Output_Filename as a tab-delimited file.

Fast_ExAC_API_Mining(numThreads): 
Handles the retrieval of allele frequency from the ExAC Project API asynchronously.
This function calls numThreads number of threads (default 200) and place them in a queue.
Each thread sends a request to the ExAC server for allele frequency.
The order of request retrieval is retained by writing to a hash table.
This creates an explicit mapping between each variant ID and allele frequency, independent of order of calculation.
If allele frequency for a variant is not available, "Not Available" is returned for that value.
    '''


if len(NumArguments) == 1:
    Error_NoArguments()
elif len(NumArguments) == 2 and NumArguments[1] == "help":
	Help_Text()
elif len(NumArguments) == 2 and NumArguments[1] != "help":
    Error_NoArguments()
elif len(NumArguments) == 3:
	import TakeHomeChallenge
	VCF = TakeHomeChallenge.VCF_Annotation()
	VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2])
elif len(NumArguments) == 4:
	import TakeHomeChallenge
	VCF = TakeHomeChallenge.VCF_Annotation()
	if NumArguments[3] == "True" or NumArguments[3] == "1":
		VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],True)
	elif NumArguments[3] == "False" or NumArguments[3] == "0":
		VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],False)
	else:
		print "Value of verbose set to default 'False'"
		VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],False)
elif len(NumArguments) == 5:
	try:
		numThreadsArgument = int(NumArguments[4])
	except ValueError:
		Error_NoArguments()
	else:
		import TakeHomeChallenge
		VCF = TakeHomeChallenge.VCF_Annotation()
		if NumArguments[3] == "True" or NumArguments[3] == "1":
			VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],True,numThreadsArgument)
		elif NumArguments[3] == "False" or NumArguments[3] == "0":
			VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],False,numThreadsArgument)
		else:
			print "Value of verbose set to default 'False'"
			VCF.Run_Full_Annotation(NumArguments[1],NumArguments[2],False, numThreadsArgument)