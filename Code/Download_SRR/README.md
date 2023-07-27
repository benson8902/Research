Some bash file need to make corrections every time!

Step 1 Download SRR file 
  1. prefetch --option-file filename.txt
  or
  2. nohup prefetch --option-file filename.txt > nohup.out 2>&1 &

Step 2 SRA -> fastq
  1. nohup fastq-dump --split-files --gzip SRR13482553

Step 3 
