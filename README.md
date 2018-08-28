# TakeHomeAssignment

This is the take home assignment for the Tempus Bioinformatics Technical Challenge. Coded in Python 2.7 using the Anaconda distribution.

TakeHomeChallenge.py is the file that contains the objects and functions associated with building the annotation.
Run_VCF_Annotation.py is a separate file that is designed to be executable by command line.

Dependencies: Uses the python packages pandas and requests. All other packages are part of the Python Standard Library.

How to use:
```python Run_VCF_Annotation.py Input_VCF_Directory Output_Directory Verbose numThreads```

Input_VCF_Directory: The directory of the VCF file to be annotated.
Output_VCF_Directory: The directory of the annotated file (outputs as tab delimited file).
Verbose (optional): Boolean, default False. When Verbose = True, all information in the original vcf file is retained in the annotation. When Verbose = False, only the chromosome, position, reference allele, alternative allele, and quality values are retained.
numThreads (optional): Int, default 200. The number of threads permissible to open when retrieving information from the ExAC Project API.

Information on program (also accessible by running ```python Run_VCF_Annotation.py help```):



## VCF_Annotation

Input: A directory to a vcf file (Input_VCF)
Output: An annotated vcf written to a specified directory (Output_Filename)

The vcf is annotated with the following information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.)  - already represented by TYPE in vcf file's INFO column
2. Depth of sequence coverage at the site of variation - already represented by DP category in vcf file's INFO column
3. Number of reads supporting the variant - already represented by AO category in vcf file's INFO column
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


## Important Functions within VCF_Annotation:

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
