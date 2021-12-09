## AS_scripts
We followed [this](https://github.com/armartin/ancestry_pipeline) pipeline by Dr. Alicia Martin for our analyses of local ancestry inference, global ancestry calculation, and ancestry-specific PCA. But, we made little modification (those with a comment "by dangliu" at the end) on the **shapeit2rfmix.py**, **lai_global.py**, and **make_aspca_inputs.py** for our analyses. I will explain the modification as following:  
*  **shapeit2rfmix.py**  
I modified line 45:  
from ```if header0 != ['ID_1', 'ID_2', 'missing', 'father', 'mother', 'sex', 'plink_pheno']:```  
to ```if header0 != ['ID_1', 'ID_2', 'missing']:```  
As the ```.sample``` output file from the version of SHAPEIT (v2) I used only have these three columns (ID_1, ID_2, missing) in the header.  
* **lai_global.py**  
I added ```pops_s = pops.split(",")``` to specifically spearate the input reference population labels by comma, which gives a tidier output format for me. So, every related line was modified accordingly, including line 24-25, line 34-35, line 44-45, and line 50.  
*  **make_aspca_inputs.py**  
The major modification is to make the ```--keep``` argument work with whoever individual I would like to keep, not just limited to reference groups. To do this, line 69 has been rivsed from ```ind_order = all_inds``` to ```ind_order = all_inds[:]``` (which makes ind_order as a copied list of all_inds) to prevent  
```
for ind in all_inds:
            if ind in keep:
                ind_order.append(ind)
```    
from keep adding ind in to the loop as the script thinks ```ind_order = all_inds```. Then, the loops at line 84-94 and line 122-131 were modified to make the ```--keep``` argument also work for admixed targets.  
One minor modification is to change line 111 ```fbk_chunked = chunker(fbk_line, 3)``` to ```fbk_chunked = chunker(fbk_line, 2)```, as I only used two reference sources.  

