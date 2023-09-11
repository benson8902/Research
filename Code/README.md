**Step 1 Preprocessing** 
>    Filter preocessed data
>    Fill in age, gender, ID into h5ad file
>    ★

**Step 2 SRA -> fastq**

  `Excute 『fastq_dump.sh』`
  
    ex:
    1. nohup fastq-dump --split-files --gzip SRR*******

**Step 3 fastqc : Quality Check**

  `Excute 『fastqc.sh』`

>  ★If there is SRR******_1.fastq.gz, SRR******_2.fastq.gz, SRR******_3.fastq.gz after fastq, you need to check which one is Read1, Read2 and Index1.
>
>  ★Even though only have SRR******_1.fastq.gz and SRR******_2.fastq.gz also need to check which one is Read1 and Read2.
> 
>  ★After check fastqc is correct, rename fastq file
>  
  `Excute 『rename_fastq.sh』`
  
      ex:
      1. fastqc -t 16 -o fastqc SRR******_1.fastq.gz
        -t => how many threads you want to use in CPU
