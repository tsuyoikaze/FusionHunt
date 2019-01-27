import os
os.system("prefetch -v [SRRNUMBER]")
os.system("fastq-dump --outdir /pylon5/mc5frap/siweixu/project --split-files ~/ncbi/public/sra/[SRRNUMBER].sra")
os.system("rm -rf ncbi")
os.system("chmod g+w *.fastq")
