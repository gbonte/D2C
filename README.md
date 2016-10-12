# D2C

Code related to the D2C algorithm introduced in the JMLR paper 
"From Dependency to Causality: A Machine Learning Approach" 


Details in http://jmlr.org/papers/v16/bontempi15a.html



To use it in your R code

- library(devtools)

- install_github("gbonte/D2C")

- require(D2C)


---------------

## Trained D2C models

For sake of space the directory /data contains D2C models trained with small number
of DAGs.

Larger D2C models are contained (in the Rdata format) in the public dropbox folder

https://www.dropbox.com/sh/t1nt236ssj7oovy/AAAVITlfRwmYVZOT6BC114-wa?dl=0

You can  load them from R with the command

load(url("https://dl.dropboxusercontent.com/u/15579987/D2Cdata/namefile.RData"))



—-----------------------
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


 
Director

Interuniversity Institute of Bioinformatics in Brussels (IB)²

Office Phone: +32-2-650 59 43

http://ibsquare.be 
