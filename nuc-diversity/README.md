
## Workflow

1) Compute frequencies using PLINK

plink --bfile ... --freqx  --out snpdata

2) Create a position file (e.g. from PLINK BIM file)
in the twwo-column format SNPID POS

awk '{print $2,$4}' ...bim > snpdata.pos

3) Run

Rscript diversity-script.R snpdata.frqx snpdata.pos 

