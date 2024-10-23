##!/usr/bin/env bash

#BSUB -q long
#BSUB -P re_gecip_renal
#BSUB -o /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/cancer_%J.stdout
#BSUB -e /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/cancer_%J.stderr
#BSUB -J cancer_controls
#BSUB -cwd /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes
#BSUB -R "rusage[mem=20000]"

#method to create mixed ancestry case controls for ancestry matched control data
#first take all proband PKDs and unrelated, non-HPO, non-HES controls and run through KING
module load bio/KING/2.2.4
module load bio/PLINK/1.9b_4.1-x86_64
module load lang/Python/3.7.4-GCCcore-8.3.0

plink -bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep-fam /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/all_cystic_platekeys_v15.txt \
--allow-no-sex --make-bed --out cystic_cases
#for stone cohort this doesn't remove anyone
king -b cystic_cases.bed --related --degree 2

#extract unrelated cases
plink --bfile cystic_cases --remove king.kin0 --allow-no-sex --make-bed --out unrelated_cystic_cases
awk '{print $1, $1}' cystic_cases.fam > unrelated_cystic_cases_ids.txt
awk '{print $1, $1}' v15_unaffected_no_renal_via_hpo_ppexplorer_controls.txt > v15_unaffected_no_renal_via_hpo_ppexplorer_controls_ids.txt

#awk '{print $1, $1}' unrelated_nonrenal_controls_platekeys.txt > tmp.txt && mv tmp.txt unrelated_nonrenal_controls_platekeys.txt
#merge with unrelated controls
cat  unrelated_cystic_cases_ids.txt v15_unaffected_no_renal_via_hpo_ppexplorer_controls_ids.txt > unrelated_cystic_rd_controls_ids.txt
#extract from plink as case/control
plink --bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep-fam unrelated_cystic_rd_controls_ids.txt \
--allow-no-sex \
--make-bed \
--out unrelated_cystic_rd_controls

#run king
king -b unrelated_cystic_rd_controls.bed --related --degree 2
#remove unrelated controls without removing cases using Catalin's script
#awk '{print $1,$1}' cancer_platekeys_v15_no_renal_bladder_etc.txt > tmp.txt && mv tmp.txt cancer_platekeys_v15_no_renal_bladder_etc.txt
/re_gecip/renal/oalavijeh/projects/phenotypes/scripts/remove-related-controls.py unrelated_cystic_cases_ids.txt v15_unaffected_no_renal_via_hpo_ppexplorer_controls_ids.txt king.kin0 to_keep.txt
#create final data set
plink --bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep to_keep.txt \
--allow-no-sex \
--make-pheno unrelated_cystic_cases_ids.txt '*' \
--make-bed \
--out /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/final_unrelated_cystic_rd_controls
#run covariate analysis
plink --bfile /re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/final_unrelated_cystic_rd_controls --pca 10 --allow-no-sex
#merge these to IDs
