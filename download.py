import os
f = open("download_srr_series.txt")
files = map(str.strip, f.readlines())
for file in files:
    os.system("prefetch -v %s", file)
    os.system("fastq-dump --outdir /pylon5/mc5frap/siweixu/project --split-files ~/ncbi/public/sra/%s", file + ".sra")
    os.system("rm -rf ~/ncbi")
os.system("chmod g+w *.fastq")
