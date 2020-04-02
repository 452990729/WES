/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/MergeData.py -i All.raw.snp.filter.hg19_multianno.txt.intervar -a ../1.Annovar/All.raw.snp.filter.anno.hg19_multianno.txt -o All.snp.interver.annovar.xls
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/AddAnno.py -i All.snp.interver.annovar.xls -o All.snp.finalanno.xls
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/FilterIntervar.py -i All.snp.finalanno.xls -o ./
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/SepSampleF.py -i FinalsnpResult.txt -v ../../4.Mutation/All.raw.s -o FinalSampleINDEL.txt
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/SepSampleF.py -i FinalindelResult.txt -v ../../4.Mutation/All.raw.indel.filter.vcf -o FinalSampleINDEL.txt
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/SepSampleF.py -i All.snp.finalanno.xls -v ../../4.Mutation/All.raw.snp.filter.vcf -o AllSampleSNP.txt
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/SepSampleF.py -i All.indel.finalanno.xls -v ../../4.Mutation/All.raw.indel.filter.vcf -o AllSampleINDEL.txt
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/IntergrateMutation.py -s AllSampleSNP.txt -i AllSampleINDEL.txt -o AllSampleMutation.txt
/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Lib/3.Intervar/IntergrateMutation.py -s FinalSampleSNP.txt -i FinalSampleINDEL.txt -o FinalSampleMutation.txt
