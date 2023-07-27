Some bash file need to make corrections every time!

Step 1 Download SRR file 

    ex:
    1. prefetch --option-file filename.txt
    or
    2. nohup prefetch --option-file "filename.txt" > nohup.out 2>&1 &

Step 2 SRA -> fastq

    ex:
    1. nohup fastq-dump --split-files --gzip SRR*******

Step 3 fastqc

  Quality check
    if there is SRR******_1.fastq.gz, SRR******_2.fastq.gz, SRR******_3.fastq.gz after fastq, you need to check which one is Read1, Read2 and Index1.
    
      ex:
      1. fastqc -t 16 -o fastqc SRR******_1.fastq.gz
        -t => how many threads you want to use in CPU
