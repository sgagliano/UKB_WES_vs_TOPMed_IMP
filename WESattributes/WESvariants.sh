#DESCRIPTION: print out WES variants

bcftools query -f "%CHROM\t%POS\t%REF","%ALT\n"  /net/hunt/disk2/UK_Biobank/Exome/ukb_spb_exm_chrall_v1.vcf.gz |bgzip -f > ../output/WES.markers.txt.gz

