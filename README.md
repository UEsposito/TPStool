# Temporal Population Structure (TPS)
Matlab source code to run the TPS tool

TPS is a method to date ancient genomes solely from their genomic sequence.
In the following we provide the code to run the TPS tool reported by Esposito et al.'s study

There are 3 folders:
-  data/ contains an excel file with the data need by TPS to run (i.e. reference individuals in the file TPS_reference_ids.xlsx) and an excel file with sample data to run TPS.
-  src/ contains the Matlab scripts necessary to run TPS.
-  results/ contains the output of TPS when run on the sample data.

To run TPS do the following:
   1. Copy the folder structure to your computer.
   2. Execute the function TPS_main.m with an input file, for example, TPS_main('../data/TPS_test_ids.xlsx'). This is the main function which loads the excel files in data/ and outputs the results to results/.
   3. The running time is <10s per sample.

Wednesday 30th October 2019
Umberto Esposito
umberto.espos@gmail.com
