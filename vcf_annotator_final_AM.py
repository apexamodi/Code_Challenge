import requests
import json
import sys

'''The goal of this script is to read in any VCF file and get the following information about the variant:

a. type of variation 
b. the effect of the variation from ExAC
c. depth of sequence coverage at the variation site
d. number of reads supporting the variant
e. percentage of reads supporting the variant versus those supporting reference reads 
f. allele frequency of variant from ExAC'''

#function to pull specified data from the INFO column of the VCF file 
def extract_variant_info(phrase,variant_info,allele_count,num):
	capture_info = list(filter(lambda x: phrase in x, variant_info))[num].split("=")[1]
	if len(capture_info) == 0:
		capture_info = "-1"
	else:
		capture_info = capture_info.split(",")[allele_count]
	return capture_info

#function to pull data from the ExAC REST API
def call_exac_api(link,variant_name_array):
	#call Exac API using post to get variant info in bulk
	response = requests.post(link, data = json.dumps(variant_name_array))
	if response.status_code != 200:
		raise ApiError('API error: {}'.format(response.status_code))
	else:
		exac_info = response.json()
	return exac_info

def main(file_name):

	vcf_file = open(file_name,"r")

	#create a dict based on the position and variant alleles: same format as ExAC input  
	vcf_info = {}

	#array to store variant names in order to do bulk Exac API request
	variant_name_array = []

	# column numbers that we will pull data from
	ref_allele = 3
	alt_allele = 4
	annotation_column = 7


	for line in vcf_file.readlines():

		line = line.strip()
		
		#skips header lines - Variants start at the line after #CHROM.
		if line[0] != "#":
			#split line into an array of strings
			split_line = line.split("\t")

			#almost all the information about the variant is in column called info 
			
			variant_info = split_line[annotation_column].split(";")

			#given that some variants have multiple alleles (esp in the case of indels), we will perform annotation for all using a loop
			allele_count = 0

			#variant type info
			var_type = extract_variant_info("TYPE=",variant_info,allele_count,0)

			#depth of coverage for the variant
			depth = extract_variant_info("DP=",variant_info,allele_count,0)
			depth = int(depth)

			#number of reads supporting the reference allele; take second term in array, picks up both PRO and RO
			ref_reads = float(extract_variant_info("RO=",variant_info,allele_count,1))
			percent_ref = round(ref_reads/depth,3)

			for allele in split_line[alt_allele].split(","):

				#name of variant in chr-pos-refAllele-varAllele format
				variant_name = "-".join(split_line[0:2])+"-"+split_line[ref_allele]+"-"+allele
				variant_name_array.append(variant_name)

				#number of reads supporting the variant allele
				variant_reads = float(extract_variant_info("AO=",variant_info,allele_count,0))

				#allele frequency based on data 
				internal_freq = float(extract_variant_info("AF=",variant_info,allele_count,0))

				percent_var = round(variant_reads/depth,3)
				
				internal_freq= round(internal_freq,3)

				vcf_info[variant_name] = [variant_name,var_type,depth,variant_reads,ref_reads,percent_var,percent_ref,internal_freq]
				allele_count+=1

	freq = call_exac_api("http://exac.hms.harvard.edu/rest/bulk/variant/variant",variant_name_array)
	effect = call_exac_api("http://exac.hms.harvard.edu/rest/bulk/variant",variant_name_array)

	for variants in variant_name_array:

		#check that allele freq information is even available for the variant from Exac
		try:
			variant_freq = round(freq[variants]["allele_freq"],3)
		except KeyError:
			variant_freq = "NA"

		vcf_info[variants].append(variant_freq)

		#Consequence of all variants is not necessarily known 
		if effect[variants]["consequence"] is None or len(effect[variants]["consequence"]) == 0:
			variant_effect = "Unknown"
		else:
			try:
				#obtain first effect associated with variant as consequences appear to be already ranked by how deleterious
				variant_effect = next(iter(effect[variants]["consequence"]))
			except KeyError:
				variant_effect = "NA"

		vcf_info[variants].append(variant_effect)

	new_file = open("vcf_annotated_final.csv","w")

	#enter the headers for the new annotated vcf file:

	header = ["Variant (chr-pos-referenceAllele-variantAllele)","VariantType","DepthOfCoverage","VariantReads","ReferenceAlleleReads","PercentageOfReads-Variant","PercentageOfReads-Ref","AlleleFreqVCF","AlleleFreqExac","VariantConsequence"]
	header = ",".join(header)+"\n"

	new_file.write(header)

	for names in vcf_info:

		string_line = [str(ints) for ints in vcf_info[names]]
		new_line = ",".join(string_line)+"\n"
		print(new_line)
		new_file.write(new_line)

	vcf_file.close()
	new_file.close()

if __name__ == "__main__":
	#Enter name of vcf file when you run the script.
	if len(sys.argv) < 2:
		print("Error: Please provide a vcf file name.")
	else:
		main(sys.argv[1]) #Challenge_data.vcf
