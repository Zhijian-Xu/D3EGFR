## Introduction
As key oncogenic drivers in non-small cell lung cancer (NSCLC), various mutations of epidermal growth factor receptor (EGFR) with variable drug sensitivities have been the major obstacle for precision medicine. For the purpose, we built a database, namely D3EGFRdb, with the clinicopathologic characteristics and drug responses of 1,339 patients harboring EGFR mutations via literature mining. Besides, we developed a deep learning-based prediction model, namely D3EGFRAI, for drug sensitivity prediction of new EGFR mutation-driven NSCLC. The D3EGFR contained functions above is freely accessible at https://www.d3pharma.com/D3EGFR/index.php.

## Suggestions
Users are recommended to search for EGFR mutation patient cases and predict drug sensitivity through the D3EGFR website. The data sets and source code required to implement the D3EGFR website can be obtained in this interface.


## Usage of the source code
### First step
Install DeepPurpose. Please refer to this link for details (https://github.com/kexinhuang12345/DeepPurpose)

### Second step
```
python drugResponse_final.py -m [mutation]
```

