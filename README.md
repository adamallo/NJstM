# NJstM
Modified version of NJst trying to better accommodate missing data.

This R script implements a modified version of the species tree reconstruction method NJst (Liu and Yu, 2011), and relies on the [phybase R package](https://faculty.franklin.uga.edu/lliu/content/phybase?).

Installation:
-------------
You only need to install R and the [phybase R package](https://faculty.franklin.uga.edu/lliu/content/phybase?) in your computer. To install R follow the common procedures of your SO. Phybase can be directly downloaded from [https://faculty.franklin.uga.edu/lliu/content/phybase?](https://faculty.franklin.uga.edu/lliu/content/phybase?) and installed from R with this command `install.packages(path_to_file, repos = NULL, type="source")`.
You can directly grab the new and modified functions from this script and use them in your own R scripts to estimate species trees, or use Rscript to use it as a sort of standalone program(see below).

Usage:
-----
`Rscript njstm.r treefile mapping method outputfile`  
* **treefile**: File with all your gene trees in phylip format.
* **mapping**: gene-copy/species mapping file.
* **method**: NJstM algorithm you want to use.
* **outputfile**: name of the reconstructed species tree (it will be saved in phylip format).

Algorithms:
----------
We implemented two different modifications of the original phybase package NJst function, called "original" and "reweighed". Moreover, we also keep the original implementation by Liu ("Liu"):

* **Original**: directly implements the algorithm as described in the paper, **which actually differs from the NJst function distributed with the phybase package when the number of individuals per species is not constant**.
 
* **Reweighed**: reweights the average gene-tree internode distances (AGID) prior to the usage of the NJ reconstruction methods in order to give the same weight to every gene tree. The original algorithm indirectly gives more weight to gene trees with more individuals per species.

* **Liu**: original function implemented by Liu in phybase, **which actually differs from the algorithm described in the NJst paper**.

**Important**: Both modified versions of the original NJst function directly discard information from missing species for certain gene families, only taking into account valid comparisons. Nevertheless, it inherits the same limitations with missing combinations across study (i.e., species A is never present when species B is present across all gene families) and the same potential problems with not random missing data. **The comparative reconstruction accuracies of these different algorithms have not been comprehensively tested yet**.

Input file format:
------------------
* **Mapping**: The gene copy/species mapping is a simple table with two columns separated by one empty space. The first column identifies the gene copy and the second the species it pertain to. Example:
```
1_A A
2_A A
1_B B
2_B B
1_C C
```
**Attention:** Due to an inherited limitation, gene-tree leaf names cannot be integers unless they represent node numbers. Species names can be either numeric or alphanumeric.

If you code the gene copies and species using the same scheme than the one in the example (genecopy=id_species), you can generate the mapping file using the following bash line on the input tree file:
`cat treefile |sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e 's/ /\'$'\n''/g' |sort|uniq|tail -n+2|gsed "s/.*\_\(.*\)$/& \1/" > mapping`
