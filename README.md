# NJstM
Modified version of NJst to better accommodate missing data.

This R script implements a modified version of the species tree reconstruction method NJst (Liu and Yu, 2011), and relies on the [phybase R package](https://faculty.franklin.uga.edu/lliu/content/phybase?).

Installation:
-------------
You only need to install R and the [phybase R package](https://faculty.franklin.uga.edu/lliu/content/phybase?) in your computer. To install R follow the common procedures of your SO. Phybase can be directly downloaded from [https://faculty.franklin.uga.edu/lliu/content/phybase?](https://faculty.franklin.uga.edu/lliu/content/phybase?) and installed from R with this command `install.packages(path_to_file, repos = NULL, type="source")`.
You can directly grab the NJstM and pair.dist.nofreq.dm functions from this script and use them in your own R scripts to estimate species trees, or use Rscript to use it as a sort of standalone program(see below).

Usage:
-----
Rscript njstm.r treefile mapping method outputfile  
*treefile: File with all your gene trees in phylip format.
*mapping: gene-copy/species mapping file.
*method: NJstM algorithm you want to use.
*outputfile: name of the reconstructed species tree (it will be saved in phylip format).

Algorithms:
----------
We implemented two different modifications of the original phybase package, called "original" and "reweighed". The first, directly implements the algorithm as described in the paper, which actually differs from the NJst function distributed with the phybase package when the number of individuals per species is not constant. The second, reweights the distances prior to the usage of the NJ reconstruction methods in order to give the same weight to every gene tree, since the original algorithm gives more weight to gene trees with more individuals per species. Both modified versions of the original NJst function can handle missing species in certain gene families much better than the original, inheriting the same limitations with missing combinations across study (i.e., species A is never present when species B is present across all gene families).

Input file format:
------------------
*Mapping: The gene copy/species mapping is a simple table with two columns separated by one empty space. The first column identifies the gene copy and the second the species it pertain to. Example:
  1_A A
  2_A A
  1_B B
  2_B B
  1_C C

