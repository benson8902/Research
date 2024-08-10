**Step 1 Preprocessing** 
>    1. Filter preocessed data
>
> 
>    2. Fill in age, gender, ID into h5ad file
>
> 
>    ★Confrim each data situation

**Step 2 Using R to Annotation cell type**
>    1. Cause raw data doesn't annotation.
>
>       
>    2. We don't need to normalize before excute Azimuth, it will normalize. 


**Step 3 Remove Batch Effect**
>   In order to decrease bias from different samples, we use bbknn to batch correction.


**Step 4 Calculate Cell Dispersion**
>   Through Euclidean calculate cell dispersion, please see the archives for details.


**Step 5 Statistical test**
>   Using gene and age to do Pearson correlation.


**Step 6 GO & KEGG** 
>    BioManager need to check whether the permission is gonna Enable system administrator rights, or you may get warning.
>
>    If you need to download GitHub package, need to get token. And the token has an expiration date.
>
>    Before conduct GO & KEGG we need to choose Pearson correltaion in Negative and Positive.


**Step 7 XGBR** 
>    We choose mean and standard of cell dispersion to build the Predicting age model.


**Note**
>    1. CITE-seq => 抗原
> 
>    2. GEX => Gene Expression
>
>    3. Seurat 可以處理CITE-seq
>
>    4. Scanpy 不可以處理CITE-seq
