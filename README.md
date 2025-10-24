Suite of Classification Algorithms based on the following article: 

Trajan Murphy*, Akshunna S. Dogra*, Hanfeng Gu, Caleb Meredith, Mark Kon, Julio En-
rique Castrillion-Candas*, for the Alzheimer’s Disease Neuroimaging Initiative, FINDER:
Feature Inference on Noisy Datasets using Eigenspace Residuals, 2025. (*Equal Contributions, correspondence should be directed to: adogra@nyu.edu).

Below is the step by step instruction to reproduce the results as shown.

Step 1. Download Data (instructions are at the bottom in the %%Download data section)

Step 2. Open MATLAB, and make sure that the current path in MATLAB is (your download path)/source. 
Type 'paths' (without quotes) in the command window to add all folder subpaths to current path

Step 3. Open InitializeParameters.m. Fill out fields in InitializeParameters. Guidelines: 

        Set parameters.data.validationType = 'Kfold' and parameters.Kfold = 1 to perform LPOCV. 
        set parameters.data.path  = (the location of your data file). If you download your datasets to the 'data' folder, you may put the empty string ''.
        set parameters.data.normalize = 1
        set parameters.multilevel.chooseTrunc = 'false';
        set parameters.parallel.on = true to run parallel loop (only if you have the parallel processing toolbox)
        set parameters.gpuarray.on = true to run on GPU (only if you have GPU), it is not recommended to set both 
                parameters.parallel.on and parameters.gpuarray.on to true

	set parameters.snapshots.k1 = (The truncation parameter for each data set). We use the following 
		ADNI: 5
		CSF: 8
		newAD: 8
		GCM: 39

	set parameters.data.label = 'Plasma_M12_ADCN' (ADNI AD vs. CN)
			 	    'Plasma_M12_CNLMCI' (ADNI CN vs. LMCI)
				    'Plasma_M12_ADLMCI' (ADNI AD vs. LMCI)
				    'SOMAscan7k_KNNimputed_AD_CN' (CSF)
				    'newAD' (newAD)
				    'GCM' = GCM;

	set parameters.multilevel.Mres = 20:20:140 (ADNI)
					 700:700:7000 (CSF and newAD)
                                         1600:1600:16000 (GCM)


Step 4. (Optional Batch Processing). 

    To perform batch processing of results, you may run CompMultiSVM2 instead of CompMultiSVM (Step 5.)
     Be sure to edit the CompMultiSVM2 file first. On line 12 of CompMultiSVM2.m, you will see the assignment
        
            D = methods.all.ValuesTable('Name', 'Value',...)
            

    The arguments come in name-value pairs. Each Name can be any character array, but each Value is a 1 x n 
    cell array containing the arguments that you want to iterate over in the batch processing. Immediately in the 
    for irow = 1:height(D) loop make the corresponding assignments to the parameters struct. An example is given in the
    current CompMultiSVM2.m file 

Step 5. 
        If you want to generate results one by one, simply type CompMultiSVM in the command window. This will generate
        the results based off of your user specified fields in the InitializeParameters.m file 

Step 6. Check the results

      The results including table of accuracy and AUC, and plots in the paper will be stored within 
      /(your download path)/results/Manual_Hyperparameter_Selection

	The results structure contains the following fields: 
	'array': 5-D array containing indices related to the (i,j)th iteration of LPOCV, the actual class value, predicted class value, and raw machine score
	'notes': string array explaining each of the dimensions in 'array'
	'DimRunTime': defunct field (ignore)
	'AUC': AUC obtained for each value of Mres for MLS and ACA, and for each learner in Benchmark
	'ROCs': ROC curves corresponding to each AUC value 
	'accuracy': same as AUC but includes the accuracy
	'accuracyBalanced':contains the balanced accuracy instead of the accuracy. This is the same as accuracy for LPOCV
	'run_time': total elapsed time to create result structure
	'creation_time': date and time at which result structure was saved.
	



%%========================
%% Downloading Data Sets
%%========================

Step 1. Get access to confidential ADNI data
      Visit the source website https://ida.loni.usc.edu/login.jsp?project=ADNI.

      Complete the data use agreement and submit your application.

      Once approved, you'll receive login credentials for the ADNI Image & Data Archive (IDA).

Step 2. Download plasma data 
      Choose 'select study' to be ADNI. In the "Search & Download" dropdown menu, select "Study Files".

      In the sidebar on the left, choose "Biospecimen" -> "Biospecimen Results".

      Find and download "Biomarkers Consortium Plasma Proteomics Project RBM Multiplex Data and Primer (Zip file)" 
      
      From the folder, extract "adni_plasma_qc_multiplex_11Nov2010.csv" and save to /data/

================== Step 3. Download phenotype data ==================================

      Also in the "Search & Download" section, find "ADNIMERGE - Packages for R" and download "ADNIMERGE_0.0.1.tar.gz".

      run the following code in R:

        install.packages("Hmisc")
        install.packages("ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
        library(ADNIMERGE)
        data("adnimerge")
        m12 <- subset(adnimerge, VISCODE=='m12')
        bl <- subset(adnimerge, VISCODE=='bl')
        write.csv(m12, "/data/adni_phenotype_m12.csv", quote = F, row.names = F)
        write.csv(bl, "/data/adni_phenotype_bl.csv", quote = F, row.names = F)
  
================== Step 4. Generate the data =====================================
      
      Now we have original data ready for use:
        "/data/adni_plasma_qc_multiplex_11Nov2010.csv"
        "/data/adni_phenotype_m12.csv"

      Go to  /source/ and type paths
      Type PrepADNI in the command window
      The binary datasets will be stored in /data/

=================================

CSF

Step 1. and Step 2 are the same as for ADNI data

Step 3. Download CSF data

      Choose 'select study' to be ADNI. In the "Search & Download" dropdown menu, select "Study Files".

      In the sidebar on the left, choose "Biospecimen" -> "Biospecimen Results".

      Find and download "CruchagaLab CSF SOMAscan7k Protein matrix postQC"
     
      Save the file "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"  to /data/

Step 4. Generate the data
     
      Now we have original data ready for use:
        "/data/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"
        "/data/adni_phenotype_bl.csv"

      Then run PrepCSF.m. The binary datasets will be stored in /source/data
===================================

newAD 
====================================
GCM 

This data set is already included in the 'data' folder; in InitializeParameters.m you may therefore set parameters.data.path = ''
====================================
Remote Sensing

Refer to Deforest_Read_Me_Data-2.pdf
======================================


