/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Software/snakemake\
    --cluster "qsub -cwd -V -b y -S /bin/bash -o {log.o} -e {log.e}"\
    -j 30 -s snakefile.txt --printshellcmds --latency-wait 10
