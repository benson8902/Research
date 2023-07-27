**Some bash file need to make corrections every time!**

**First of all:**

    Need to check whether sratoolkit working.
    『export PATH=$PATH:$HOME/sratoolkit.3.0.1-ubuntu64/bin』

**Step 1 Download SRR file** 

    ex:
    1. prefetch --option-file filename.txt
    or
    2. nohup prefetch --option-file "filename.txt" > nohup.out 2>&1 &

**Step 2 SRA -> fastq**

    ex:
    1. nohup fastq-dump --split-files --gzip SRR*******

**Step 3 fastqc**

  Quality check

  ★If there is SRR******_1.fastq.gz, SRR******_2.fastq.gz, SRR******_3.fastq.gz after fastq, you need to check which one is Read1, Read2 and Index1.
    
      ex:
      1. fastqc -t 16 -o fastqc SRR******_1.fastq.gz
        -t => how many threads you want to use in CPU

**Step 4 Do 10X Genomics Cell Ranger**

    following:
    https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

  ★If there are multiple samples of a single individual then we excute cellranger aggr

**Notice:**

    SRR_Acc_List.txt need to add new line, or the last file will not be excute.
