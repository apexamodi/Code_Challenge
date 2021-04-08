# Code_Challenge

The script tempus_vcf_annotator_final.py can be used to annotate a given VCF file with the following information:

* type of variation
* the predicted consequence of the variation from ExAC
* depth of sequence coverage at the variation site
* number of reads supporting the variant
* percentage of reads supporting the variant versus those supporting reference reads 
* allele frequency of variant from ExAC

The script requires you to have python installed and can be run as follows in your terminal. You must provide a vcf file as input as shown: 

python tempus_vcf_annotator_test.py Challenge_data.vcf

Required python modules: 
1) sys
2) requests
3) json

This script makes requests from the ExAC REST API and will require access to the internet. 

Output: The script will output an annotated csv file called: vcf_annotated_final.csv The output file will contain 10 columns, the header explains the information contained in each column. 
