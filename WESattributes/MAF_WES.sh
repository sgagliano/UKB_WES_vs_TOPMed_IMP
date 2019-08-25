#DESCRIPTION: compute minor allele frequency for WES variants (only include the samples who also were TOPMed imputed

dir=/net/hunt/disk2/UK_Biobank/Exome

plink-1.9 --bed $dir/ukb_evc_chr1_v1.bed --bim $dir/ukb_spb_exm_chrall_v1.bim --fam $dir/ukb24460_evc_chr1_v1_s49959.fam --keep /net/inpsyght/disk2/sarahgag/UKB500/VEP-TOPMedIMPUTATION/VEP-TOPMed/EXOME/IDs_TOPMedImputed_EXOME_overlap-4PLINK --freq --out ../output/ukb_spb_exm
