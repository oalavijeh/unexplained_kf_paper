# script to pull out APOL1 variants from cases and control
#G1 is defined as rs73885319(G) and rs60910145(G)
#G2 is defined as rs71785313(delTTATAA)

# chr22:36265860 A>G this
# chr22:36265988 T>G and this together are G1
# chr22:36265996-36266005 delTTATAA this is G2

module load bio/BCFtools/1.11-GCC-8.3.0

#bcftools merge -r chr22:36265860,chr22:36265988,chr22:36265996-36266005 -l all_african_unexplained_cases_gvcf_paths.txt -m none -Oz -o apol1_africans.vcf.gz
#expanded is how you get G2
bcftools merge -r chr22:36265860,chr22:36265988,chr22:36265990-36266010 -l  all_african_unexplained_cases_gvcf_paths.txt -m none -Oz -o apol1_africans_expanded.vcf.gz
bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%ID\t%SAMPLE:%GT\n]' -i "GT='alt'" apol1_africans_expanded.vcf.gz > esrf_apol1_african.txt
#expanded non-african cases
bcftools merge -r chr22:36265860,chr22:36265988,chr22:36265990-36266010 -l non_african_unexplained_cases_gvcf_paths.txt -m none -Oz -o apol1_non_africans_expanded.vcf.gz
bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%ID\t%SAMPLE:%GT\n]' -i "GT='alt'" apol1_non_africans_expanded.vcf.gz > apol1_esrf_non_africans_alt_expanded.txt
#controls
bcftools merge -r chr22:36265860,chr22:36265988,chr22:36265990-36266010 -l control_african_gvcf_paths.txt -m none -Oz -o control_apol1_africans.vcf.gz
bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%ID\t%SAMPLE:%GT\n]' -i "GT='alt'" control_apol1_africans.vcf.gz >  control_apol1_africans_expanded.txt 

#search a single case

#bcftools merge -r chr22:36265860,chr22:36265988,chr22:36265996-36266005 -l all_african_unexplained_cases_gvcf_paths.txt -m none -Oz -o apol1_africans.vcf.gz
#expanded is how you get G2
bcftools merge -r chr22:36265628 -l  all_african_unexplained_cases_gvcf_paths.txt -m none -Oz -o apol1_protective_african_cases.vcf.gz
bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%ID\t%SAMPLE:%GT\n]' -i "GT='alt'" apol1_protective_african_cases.vcf.gz > esrf_apol1_protective_african.txt
