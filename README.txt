Example Run:



fragment_to_infinity:
python3 fragment_to_infinity.py ebola_uniprot-proteome%3AUP000140031.fasta -s 9 -a HLA_Class_1/A.xlsx -b HLA_Class_1/B.xlsx -c HLA_Class_1/C.xlsx




mhc_analysis:
scatter -
python3 mhc_analysis.py predictions.csv -a HLA_Class_1/A.xlsx -b HLA_Class_1/B.xlsx -c HLA_Class_1/C.xlsx --scatterF

box -
python3 mhc_analysis.py predictions.csv -a HLA_Class_1/A.xlsx -b HLA_Class_1/B.xlsx -c HLA_Class_1/C.xlsx --boxF

If non-gui interface, "-s"
python3 mhc_analysis.py predictions.csv -a HLA_Class_1/A.xlsx -b HLA_Class_1/B.xlsx -c HLA_Class_1/C.xlsx --scatterF -s results.png




Addtional Information:
Due to size constraints of prediction files, I was not able to include them in the GitHub repo. If needed, I can point to the location on the GPU server where I have them saved.

To change from Broad_race to Race - use command line option "--code" while submitting to mhc_analysis.
