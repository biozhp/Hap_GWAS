gcta64 --bfile input --autosome --maf 0.05 --make-grm --out output_grm
gcta64 --reml --grm output_grm --pheno phenotype.phen --out output
