# D2C

Code related to the D2C algorithm introduced in the JMLR paper 
"From Dependency to Causality: A Machine Learning Approach" 


Details in http://jmlr.org/papers/v16/bontempi15a.html



To use it in your R code

* library(devtools); install_github("gbonte/D2C");  require(D2C)

* if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
* BiocManager::install("graph")
* BiocManager::install("Rgraphviz")
* BiocManager::install("RBGL")


---------------

## Trained D2C models

For sake of space the directory /data contains D2C models trained with small number
of DAGs.

To change the directory to the one containing scripts 

- setwd(find.package("D2C"))

To go to the directory of the  scripts 
- setwd(paste(find.package("D2C"),"scripts",sep="/"))

To list the D2C demos:
- dir(paste(find.package("D2C"),"scripts/",sep="/"),patt="demo*")




### Author 

Pr. Gianluca Bontempi

co-Head of the Machine Learning Group

Département d'Informatique

Université Libre de Bruxelles

Boulevard du Triomphe - CP212

1050 Bruxelles, Belgium

email: gbonte@ulb.ac.be

Office Phone: +32-2-650 55 91

Fax: +32 2 650.56.09

http://mlg.ulb.ac.be


