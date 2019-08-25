#split correlation output into chromsomes
gzip -dc ../output/Autosomes.txt.gz | awk '$1==1'> ../output/Autosomes-chr1.txt
cat ../output/Header ../output/Autosomes-chr1.txt | gzip -c > ../output/Autosomes-chr1.txt.gz
rm ../output/Autosomes-chr1.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==2'> ../output/Autosomes-chr2.txt
cat ../output/Header ../output/Autosomes-chr2.txt | gzip -c > ../output/Autosomes-chr2.txt.gz
rm ../output/Autosomes-chr2.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==3'> ../output/Autosomes-chr3.txt
cat ../output/Header ../output/Autosomes-chr3.txt | gzip -c > ../output/Autosomes-chr3.txt.gz
rm ../output/Autosomes-chr3.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==4'> ../output/Autosomes-chr4.txt
cat ../output/Header ../output/Autosomes-chr4.txt | gzip -c > ../output/Autosomes-chr4.txt.gz
rm ../output/Autosomes-chr4.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==5'> ../output/Autosomes-chr5.txt
cat ../output/Header ../output/Autosomes-chr5.txt | gzip -c > ../output/Autosomes-chr5.txt.gz
rm ../output/Autosomes-chr5.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==6'> ../output/Autosomes-chr6.txt
cat ../output/Header ../output/Autosomes-chr6.txt | gzip -c > ../output/Autosomes-chr6.txt.gz
rm ../output/Autosomes-chr6.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==7'> ../output/Autosomes-chr7.txt
cat ../output/Header ../output/Autosomes-chr7.txt | gzip -c > ../output/Autosomes-chr7.txt.gz
rm ../output/Autosomes-chr7.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==8'> ../output/Autosomes-chr8.txt
cat ../output/Header ../output/Autosomes-chr8.txt | gzip -c > ../output/Autosomes-chr8.txt.gz
rm ../output/Autosomes-chr8.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==9'> ../output/Autosomes-chr9.txt
cat ../output/Header ../output/Autosomes-chr9.txt | gzip -c > ../output/Autosomes-chr9.txt.gz
rm ../output/Autosomes-chr9.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==10'> ../output/Autosomes-chr10.txt
cat ../output/Header ../output/Autosomes-chr10.txt | gzip -c > ../output/Autosomes-chr10.txt.gz
rm ../output/Autosomes-chr10.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==11'> ../output/Autosomes-chr11.txt
cat ../output/Header ../output/Autosomes-chr11.txt | gzip -c > ../output/Autosomes-chr11.txt.gz
rm ../output/Autosomes-chr11.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==12'> ../output/Autosomes-chr12.txt
cat ../output/Header ../output/Autosomes-chr12.txt | gzip -c > ../output/Autosomes-chr12.txt.gz
rm ../output/Autosomes-chr12.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==13'> ../output/Autosomes-chr13.txt
cat ../output/Header ../output/Autosomes-chr13.txt | gzip -c > ../output/Autosomes-chr13.txt.gz
rm ../output/Autosomes-chr13.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==14'> ../output/Autosomes-chr14.txt
cat ../output/Header ../output/Autosomes-chr14.txt | gzip -c > ../output/Autosomes-chr14.txt.gz
rm ../output/Autosomes-chr14.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==15'> ../output/Autosomes-chr15.txt
cat ../output/Header ../output/Autosomes-chr15.txt | gzip -c > ../output/Autosomes-chr15.txt.gz
rm ../output/Autosomes-chr15.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==16'> ../output/Autosomes-chr16.txt
cat ../output/Header ../output/Autosomes-chr16.txt | gzip -c > ../output/Autosomes-chr16.txt.gz
rm ../output/Autosomes-chr16.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==17'> ../output/Autosomes-chr17.txt
cat ../output/Header ../output/Autosomes-chr17.txt | gzip -c > ../output/Autosomes-chr17.txt.gz
rm ../output/Autosomes-chr17.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==18'> ../output/Autosomes-chr18.txt
cat ../output/Header ../output/Autosomes-chr18.txt | gzip -c > ../output/Autosomes-chr18.txt.gz
rm ../output/Autosomes-chr18.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==19'> ../output/Autosomes-chr19.txt
cat ../output/Header ../output/Autosomes-chr19.txt | gzip -c > ../output/Autosomes-chr19.txt.gz
rm ../output/Autosomes-chr19.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==20'> ../output/Autosomes-chr20.txt
cat ../output/Header ../output/Autosomes-chr20.txt | gzip -c > ../output/Autosomes-chr20.txt.gz
rm ../output/Autosomes-chr20.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==21'> ../output/Autosomes-chr21.txt
cat ../output/Header ../output/Autosomes-chr21.txt | gzip -c > ../output/Autosomes-chr21.txt.gz
rm ../output/Autosomes-chr21.txt

gzip -dc ../output/Autosomes.txt.gz | awk '$1==22'> ../output/Autosomes-chr22.txt
cat ../output/Header ../output/Autosomes-chr22.txt | gzip -c > ../output/Autosomes-chr22.txt.gz
rm ../output/Autosomes-chr22.txt

#run annotate.py
for CHR in `seq 1 22`
do
python annotate.py -i ../output/Autosomes-chr${CHR}.txt.gz -c /net/topmed/working/dtaliun/TOPMed_paper_65k/VEP_LoF.02/chr${CHR}.vep.vcf.gz -o ../output/annotated-chr${CHR}.txt.gz
done
