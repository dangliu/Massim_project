## AS_scripts
We followed [this](https://github.com/armartin/ancestry_pipeline) pipeline by Dr. Alicia Martin for our analyses of local ancestry inference, global ancestry calculation, and ancestry-specific PCA. But, we made little modification on the **shapeit2rfmix.py**, **lai_global.py**, and **make_aspca_inputs.py** for our analyses. I will explain the modification as following:  
*  shapeit2rfmix.py
I modified line 45:
from
```if header0 != ['ID_1', 'ID_2', 'missing', 'father', 'mother', 'sex', 'plink_pheno']:```
to
```if header0 != ['ID_1', 'ID_2', 'missing']:```
As the ```.sample``` output file from the version of SHAPEIT (v2) I used only have these three columns (ID_1, ID_2, missing) in the header.
