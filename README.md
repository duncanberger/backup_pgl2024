<h1 align="center">
  Pneumococcal Genome Library cgMLST Typing Scheme
</h1>

![](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![](https://img.shields.io/badge/RStudio-75AADB?style=for-the-badge&logo=RStudio&logoColor=white)
![](https://img.shields.io/badge/Shell_Script-121011?style=for-the-badge&logo=gnu-bash&logoColor=white)
![](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)  
![](https://img.shields.io/badge/Repository_created:-29_May_2024-green)
![](https://img.shields.io/badge/Last_update:-29_May_2024-black)
![](https://img.shields.io/badge/PRIVATE-red)

<p align="center">
  <img width="300" height="300" src=pgl.png>
</p>

# Overview
This repository contains the code used to run analyses and generate figures for the following manuscript:  
<p align="center">
<a href="www.page.com">*Microbial Genomics*</a>
</p>  

>[!WARNING]
>Link not available yet - will be updated once manuscript goes live...  

> [!NOTE]  
> This repository contains modified code from the original repository that accompanied the pre-print publication in [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.19.571883v1). The original repository and code was authored by [Duncan Burger](https://github.com/duncanberger). To view the original repository, please click on the following links:
> |[Brueggemann Lab (forked) Repository](https://github.com/brueggemann-lab/PGL_cgMLST) | [Duncan Burger's Repository](https://github.com/duncanberger/PGL_cgMLST)|
> |-------------------------------------------------------------------------------------|-------------------------------------------------------------------------|   

# Codebook
This repository contains two main folders: `Analysis` and `Figures`. 
```mermaid
graph TD;
repo[This Repository]
dist[Distance Matrix]
tree[Phylogenetic Tree]
clust[GPSC Clustering]
clust2[Mandrake Clustering]
repo-->Analysis;
Analysis-->dist;
Analysis-->tree;
Analysis-->clust;
Analysis-->clust2;
dist-->Manuscript;
tree-->Manuscript;
clust-->Manuscript;
clust2-->Manuscript;
repo-->Figures;
Figures-->Manuscript;
dist-->Figures;
tree-->Figures;
clust-->Figures;
clust2-->Figures;
```
<details>
<summary>View folder contents</summary>
<ol>
  <li>Analysis - contains the code used to generate:</li>
  <ol>
      <li>the distance matrix</li>
      <li>the phylogenetic tree</li>
    <li>GPSC and Mandrake clustering</li>
    </ol>
  <li>Figures - contains the R code in markdown format used to generate main and supplementary figures</li>
</ol>
</details>

# License
Distributed under the GNU General Public License v3.0. Please see `LICENSE` for more information.
## Publication History
|**Publication**|**DOI**|
|-------------------------------|------|
|[*bioRxiv 2023*](https://www.biorxiv.org/content/10.1101/2023.12.19.571883v1)|![doi](https://img.shields.io/badge/DOI-https://doi.org/10.1101/2023.12.19.571883-blue)|  

>[!TIP]
>**Reference**: Jansen van Rensburg MJ, Berger DJ, Fohrmann A, Bray JE, Jolley KA, Maiden MC, Brueggemann AB. Development of the Pneumococcal Genome Library, a core genome multilocus sequence typing scheme, and a taxonomic life identification number barcoding system to investigate and define pneumococcal population structure. bioRxiv. 2023:2023-12.
# Contact
If you have any queries, suggestions or concerns, please contact [Professor Angela Brueggemann](mailto:angela.brueggemann@ndph.ox.ac.uk)  

|Repository `PRIVATE` status since: 29/05/2024|
|--------------------------------------------|
