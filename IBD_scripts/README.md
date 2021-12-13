## IBD_scripts
This directory contains R scripts for analyzing between-group IBD sharing among Oceanians (**IBD_network_bar.R**) and IBD similarity among Kula vs. non-Kula groups through space and time (**Geo_vs_IBD.R** and **Relative_IBD_similarity_stats.R**).  
***
### Library requirements
* tidyverse
* ggmap
* maps
* ggrepel
* reshape2
* geosphere
Can run this command in R to install all the packages: ```install.packages(c("tidyverse", "ggmap", "maps", "ggrepel", "reshape2", "geosphere"))```  
***
The IBD result file is **all.lPSC.Merged** which contains inferred and merged IBD between all the individuals in our dataset by [refinedIBD](https://faculty.washington.edu/browning/refined-ibd.html) and their merge-ibd-segments.jar java script. The .ibd and .hbd output files were merged all together as long Pairwise Shared Coalescent (lPSC) according to [previous study](https://github.com/halasadi/MAPS).  
And an meta information file **Kula.meta.info**, containing individual ID, group label, geo-coordinates, QC results...etc., is used to filter and label the individuals with their groups.  
The results of **IBD_network_bar.R** were used for the Fig. 7A and fig. S13 in our manuscript. The result of **Geo_vs_IBD.R** was used for fig. S16, and the results of **Relative_IBD_similarity_stats.R** were for Fig. 8 and fig. S17.