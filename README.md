## Introduction
As key oncogenic drivers in non-small cell lung cancer (NSCLC), various mutations of epidermal growth factor receptor (EGFR) with variable drug sensitivities have been the major obstacle for precision medicine. For the purpose, we built a database, namely D3EGFRdb, with the clinicopathologic characteristics and drug responses of 1,339 patients harboring EGFR mutations via literature mining. Besides, we developed a deep learning-based prediction model, namely D3EGFRAI, for drug sensitivity prediction of new EGFR mutation-driven NSCLC. The D3EGFR contained functions above is freely accessible at https://www.d3pharma.com/D3EGFR/index.php.

The D3EGFRdb database file is *D3EGFR-database.csv*, which can be used for patient case retrieval.  
The running script of D3EGFRAI is *drugResponse_final.py* for drug response prediction.

## Suggestions
Users are recommended to search for EGFR mutation patient cases and predict drug sensitivity through the D3EGFR website. The data sets and source code required to implement the D3EGFR website can be obtained in this interface.

## Usage of the source code
### First step: Install DeepPurpose. 
The drug and protein encoders were provided by *DeepPurpose*.  
Therefore, before running D3EGFR, please make sure *DeepPurpose* is installed correctly.  
Please refer to this link for how to install *DeepPurpose*: https://github.com/kexinhuang12345/DeepPurpose  
DeepPurpose doi: 10.1093/bioinformatics/btaa1005

### Second step: run D3EGFR
```
#mutation examples:
#Point mutation: L858R
#Deletion mutation: E746_A750del, V834del
#Insert mutation: D770insSVD
#Duplicate mutation: A767dupASV
#Deletion-insertion mutation: L747_P753delinsS
#Complex mutation: E746_A750del+L858R

python drugResponse_final.py -m [mutation]
```
## Other
The patient dataset and external dataset is also released in zenodo (https://zenodo.org/records/10612678). 

If you have any questions about the installation and use of D3EGFR, please contact us (zjxu@simm.ac.cn).


