## Introduction
As key oncogenic drivers in non-small cell lung cancer (NSCLC), various mutations of epidermal growth factor receptor (EGFR) with variable drug sensitivities have been the major obstacle for precision medicine. For the purpose, we built a database, namely **D3EGFRdb**, with the clinicopathologic characteristics and drug responses of 1,339 patients harboring EGFR mutations via literature mining. Besides, we developed a deep learning-based prediction model, namely **D3EGFRAI**, for drug sensitivity prediction of new EGFR mutation-driven NSCLC. The D3EGFR contained functions above is freely accessible at **https://www.d3pharma.com/D3EGFR/index.php**.

The D3EGFRdb database file is **D3EGFR-database.csv**, which can be used for patient case retrieval.  
The running script of D3EGFRAI is **drugResponse_final.py** for drug response prediction.

## Suggestions
Users are recommended to search for EGFR mutation patient cases and predict drug sensitivity through the D3EGFR website. The data sets and source code required to implement the D3EGFR website can be obtained in this interface.

## Documentation
The D3EGFR website manual can be found in https://www.d3pharma.com/D3EGFR/help.php. For local installation and use of D3EGFR, please refer to the following steps.
### Installation
#### 1.Install DeepPurpose  
> The drug and protein encoders are provided by **DeepPurpose**.  
Therefore, before running D3EGFR, please make sure **DeepPurpose** is installed correctly.  
Please refer to this link for how to install **DeepPurpose**: https://github.com/kexinhuang12345/DeepPurpose  
DeepPurpose **doi**: 10.1093/bioinformatics/btaa1005  
To install locally, we recommend to install **DeepPurpose** from **pip**:
```
conda create -n DeepPurpose python=3.6
conda activate DeepPurpose
conda install -c conda-forge notebook
pip install git+https://github.com/bp-kelley/descriptastorus 
pip install DeepPurpose
```

#### 2.Download D3EGFR code
```
git clone https://github.com/Zhijian-Xu/D3EGFR.git
```
Or you can download and unzip the D3EGFR zip file.

### USAGE
> Run the following command under Linux to predict drug response.
```
python drugResponse_final.py -m [mutation]

#[mutation] examples:
#Point mutation: L858R
#Deletion mutation: E746_A750del, V834del
#Insert mutation: D770insSVD
#Duplicate mutation: A767dupASV
#Deletion-insertion mutation: L747_P753delinsS
#Complex mutation: E746_A750del+L858R
```
## License
The code in this package is licensed under the MIT License.

## Other
The patient dataset and external dataset are also released in **zenodo** (https://zenodo.org/records/10613332). 

If you have any questions about the installation and usage of D3EGFR, please contact us (zjxu@simm.ac.cn; xbzhang@simm.ac.cn).
