Suite of Classification Algorithms based on the following article: 

Trajan Murphy*, Akshunna S. Dogra*, Hanfeng Gu, Caleb Meredith, Mark Kon, Julio En-
rique Castrillion-Candas*, for the Alzheimer’s Disease Neuroimaging Initiative, FINDER:
Feature Inference on Noisy Datasets using Eigenspace Residuals, 2025. 

(*Equal Contributions, correspondence should be directed to: adogra@nyu.edu).

Below is the step by step instruction to reproduce the results as shown.

Step 1. Download Data (instructions are at the bottom in the %%Download data section)

Step 2. Type 'paths' (without quotes) in the command window to add all folder subpaths to current path

Step 3. Fill out fields in InitializeParameters. Guidelines: 

        Set parameters.data.validationType = 'Kfold' and parameters.Kfold = 1 to perform LPOCV. 
        set parameters.data.path  = (the location of your data file) 
        set parameters.data.normalize = 1
        set parameters.multilevel.chooseTrunc = 'false';
        set parameters.parallel.on = true to run parallel loop (only if you have the parallel processing toolbox)
        set parameters.gpuarray.on = true to run on GPU (only if you have GPU), it is not recommended to set both 
                parameters.parallel.on and parameters.gpuarray.on to true

Step 4. (Optional Batch Processing). 

    To perform batch processing of results, you may run CompMultiSVM2 instead of CompMultiSVM (Step 5.)
     Be sure to edit the CompMultiSVM2 file first. On line 12 of CompMultiSVM2.m, you will see the assignment
        
            D = methods.all.ValuesTable('Name', 'Value',...)
            

    The arguments come in name-value pairs. Each Name can be any character array, but each Value is a 1 x n 
    cell array containing the arguments that you want to iterate over in the batch processing. Immediately in the 
    for irow = 1:height(D) loop make the corresponding assingments to the parameters struct. An example is given in the
    current CompMultiSVM2.m file 

Step 5. 
        If you want to generate results one by one, simply type CompMultiSVM in the command window. This will generate
        the results based off of your user specified fields in the InitializeParameters.m file 

Step 6. Check the results

      The results including table of accuracy and AUC, and plots in the paper will be stored in 
      /yourpath/TensorStochasticMachineLearning/results/



%%========================
%% Downloading Data Sets
%%========================

ADNI
========================

Step 1. Get access to confidential ADNI data

      Visit the source website https://ida.loni.usc.edu/login.jsp?project=ADNI.

      Complete the data use agreement and submit your application.

      Once approved, you'll receive login credentials for the ADNI Image & Data Archive (IDA).

Step 2. Download plasma data

      In the "Search & Download" dropdown menu, select "Study Files".

      In the sidebar on the left, choose "Biospecimen" -> "Biospecimen Results".

      Find and download "Biomarkers Consortium Plasma Proteomics Project RBM Multiplex Data and Primer (Zip file)" 
      
      From the folder, extract "adni_plasma_qc_multiplex_11Nov2010.csv" and save to /yourfolder

Step 3. Download phenotype data

      Also in the "Search & Download" section, find "ADNIMERGE - Packages for R" and download "ADNIMERGE_0.0.1.tar.gz".

      run the following code in R:

        install.packages("Hmisc")
        install.packages("/your/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
        library(ADNIMERGE)
        data("adnimerge")
        m12 <- subset(adnimerge, VISCODE=='m12')
        bl <- subset(adnimerge, VISCODE=='bl')
        write.csv(m12, "/yourfolder/adni_phenotype_m12.csv", quote = F, row.names = F)
        write.csv(bl, "/yourfolder/adni_phenotype_bl.csv", quote = F, row.names = F)

Step 4. Generate the data
      
      Now we have original data ready for use:
        "/yourfolder/adni_plasma_qc_multiplex_11Nov2010.csv"
        "/yourfolder/adni_phenotype_m12.csv"

      Open /yourpath/TensorStochasticMachineLearning/source/Modules/PrepADNI.m, 
      and change file paths for plasma and phenotype to your data location. 
      Then run PrepADNI.m. The binary datasets will be stored in /source/ADNI_data
=================================

CSF

Step 1. and Step 2 are the same as for ADNI data

Step 3. Download CSF data

      Choose 'select study' to be ADNI. In the "Search & Download" dropdown menu, select "Study Files".

      In the sidebar on the left, choose "Biospecimen" -> "Biospecimen Results".

      Find and download "CruchagaLab CSF SOMAscan7k Protein matrix postQC"
     
      Save the file "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"  to /yourfolder

Step 4. Generate the data
     
      Now we have original data ready for use:
        "/yourfolder/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"
        "/yourfolder/adni_phenotype_bl.csv"

      Open /yourpath/TensorStochasticMachineLearning/source/Modules/PrepCSF.m,
      and change file paths for plasma and phenotype to your data location.
      Then run PrepCSF.m. The binary datasets will be stored in /source/data/CSF_data
===================================

newAD 
====================================
GCM 

This data set is already included in the 'data' folder; in InitializeParameters.m you may therefore set parameters.data.path = ''
====================================
Remote Sensing

Refer to Deforest_Read_Me_Data-2.pdf
======================================

