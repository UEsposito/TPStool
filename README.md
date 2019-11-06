# Temporal Population Structure (TPS)
This is the Matlab source code necessary to run the TPS tool as a standalone tool when you wish to date a new sample.

TPS is a method that dates ancient genomes solely from their genomic sequence.
In the following we provide the code to run the TPS tool reported by Esposito et al.'s study

There are 3 folders:
-  data/ contains two excel files with the data needed by TPS to run. First, TPS_reference_ids.xlsx, this file contains the temporal components and dates for the  reference individuals. This file should not be modified for Eurasian samples, unless you want to add new references. The second excel file has the sample's temporal components which TPS analyses to infer their date. To get these temproal components, run supervised ADMIXTURE for your sample agaisnt the eight temporal populations provided in the data folder of the TPSpaper repository.
-  src/ contains the Matlab scripts necessary to run TPS.
-  results/ contains the output of TPS when run on the sample data.

To run TPS do the following:
   1. Copy the folder structure to your computer.
   2. Execute the function TPS_main.m with an input excel file, for example, TPS_main('../data/TPS_test_ids.xlsx'). This is the main function which loads the excel files in data/ and outputs the results to results/.
   3. The running time is <10s per sample.

Wednesday 30th October 2019
Umberto Esposito
umberto.espos@gmail.com
