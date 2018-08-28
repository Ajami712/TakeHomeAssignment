# coding: utf-8

from copy import deepcopy
import requests
import json
from threading import Thread
from Queue import Queue
import pandas as pd


class VCF_Annotation():
    '''
    #=======================================
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
    
    def __init__(self):
        self.Indices = []
        self.TextBody = []
        self.ColNames = []
        self.Metadata = ''
        
    def Fast_ExAC_API_Mining(self, numThreads = 200):
        
        if(isinstance(numThreads,(int, long)) == False):
            raise TypeError("Number of threads must be an integer.")
            
        URLs = ["http://exac.hms.harvard.edu/rest/variant/variant/" + i for i in self.Indices]
        concurrent = numThreads
        Map_To_Indices = {}

        def mining():
            while True:
                url = q.get()
                status, url = attemptRetrieval(url)
                mapAlleleFreqToIndex(status, url)
                q.task_done()

        def attemptRetrieval(ourl):
            try:
                response = requests.get(ourl)
                if response.status_code != 200:
                    return "error", ourl
                elif 'allele_freq' not in response.json():
                    return "Not Available", ourl
                else:
                    return str(response.json()['allele_freq']), ourl
            except:
                return "error", ourl

        def mapAlleleFreqToIndex(status, url):
            Map_To_Indices[ url.split("variant/")[-1] ] = status

        
        q = Queue(concurrent * 2)
        for i in range(concurrent):
            t = Thread(target=mining)
            t.daemon = True
            t.start()
        try:
            for url in URLs:
                q.put(url.strip())
            q.join()
        except KeyboardInterrupt:
            sys.exit(1)
            
        return Map_To_Indices
        
    def Get_VCF_Lines(self,Input_VCF):
        with open(Input_VCF,"r") as vcf_original:
            return sum(1 for line in vcf_original)
            
    def Parse_VCF(self, Input_VCF):
        
        num_lines = self.Get_VCF_Lines(Input_VCF)
        with open(Input_VCF,"r") as vcf_original:
            for i in range(0, num_lines):
                Active_Line = vcf_original.readline()
                if len(Active_Line) == 0: # End of file
                    break
                elif Active_Line[0:2] == "##": # Metadata of vcf file
                    self.Metadata += Active_Line
                elif Active_Line[0] == "#": # Header of annotation
                    Header = Active_Line.strip().split("\t")
                    Header = "\t".join(Header) + ("\tType of Variation" + 
                                                  "\tSequence Depth Coverage at Variation Site" + 
                                                  "\tNumber Of Reads Supporting Variant" + 
                                                  "\t% Reads Supporting Variant")
                    self.ColNames = Header.split("\t")
                    
                else: # Variants
                    Active_Line_Decomposed = Active_Line.strip().split("\t")
                
                    # Check how many variants are listed on the line
                    NumberOfVariants = len(Active_Line_Decomposed[4].split(","))
            
                    for j in range(0,NumberOfVariants):
                        # Presumably, we only want to retain the lines to the left of the "INFO" column from the original VCF.
                        Retained_SNP_Information = deepcopy(Active_Line_Decomposed)
                        Retained_SNP_Information[4] = Active_Line_Decomposed[4].split(",")[j]
                
                        # The rest API requires chromosome, position, reference, and alternative allele
                        Information_For_Rest_API = (Retained_SNP_Information[0] + "-" + 
                                                    Retained_SNP_Information[1] + "-" +
                                                    Retained_SNP_Information[3] + "-" +
                                                    Retained_SNP_Information[4])
                        self.Indices += [Information_For_Rest_API]
                
                    
                        # We extract sequence depth, number of reads supporting variant, and type of variation
                        INFO = Active_Line_Decomposed[7].split(";")
                        Information_We_Want = ["TYPE","DP","AO"]
                        for info in Information_We_Want:
                            Relevant_Info = [k.split("=")[1] for k in INFO if k[0:len(info)+1]==info+"="]
                            if info == "AO" and NumberOfVariants > 1:
                                Relevant_Info = [Relevant_Info[0].split(",")[j]]
                            Retained_SNP_Information += Relevant_Info
                        Retained_SNP_Information += [str(100.0 * float(Retained_SNP_Information[-1]) / float(Retained_SNP_Information[-2]))]
                        self.TextBody += [Retained_SNP_Information]
                        
    def Finalize_VCF(self, Output_Filename, Verbose = False, numThreads = 200):
        
        if(isinstance(numThreads,(int, long)) == False):
            raise TypeError("Number of threads must be an integer.")
        if(isinstance(Verbose,(bool)) == False):
            raise TypeError("Verbose variable must be a boolean.")
        
        Index_Map = self.Fast_ExAC_API_Mining(numThreads)
        
        DF = pd.DataFrame(self.TextBody,index=self.Indices,columns=self.ColNames)
        DF["Allele Frequency of Variant"] = pd.Series(Index_Map)
        
        if Verbose == False:
            try:
                DF = DF.drop('FILTER', 1)
                DF = DF.drop('INFO', 1)
                DF = DF.drop('ID', 1)
                DF = DF.drop('FORMAT', 1)
                DF = DF.drop('normal', 1)
                DF = DF.drop('vaf5',1)
            except:
                pass
            DF.to_csv(Output_Filename, header=True, index=False, sep='\t', mode='w')
        if Verbose == True:
            with open(Output_Filename,"w") as outfile:
                outfile.write(self.Metadata)
            DF.to_csv(Output_Filename, header=True, index=False, sep='\t', mode='a')
        
    def Run_Full_Annotation(self, Input_VCF, Output_Filename, Verbose = False, numThreads = 200):
        print "Parsing VCF..."
        self.Parse_VCF(Input_VCF)
        print "Fetching Data from API and finalizing VCF..."
        self.Finalize_VCF(Output_Filename, Verbose, numThreads)
        print "Done!"

