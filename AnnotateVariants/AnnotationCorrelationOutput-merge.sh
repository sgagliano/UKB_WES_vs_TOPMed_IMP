gzip -dc ../output/annotated-chr1.txt.gz > ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr2.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr3.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr4.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr5.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr6.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr7.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr8.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr9.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr10.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr11.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr12.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr13.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr14.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr15.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr16.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr17.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr18.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr19.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr20.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr21.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip -dc ../output/annotated-chr22.txt.gz | grep -v CHROM >> ../output/annotated-Autosomes.txt
gzip ../output/annotated-Autosomes.txt
