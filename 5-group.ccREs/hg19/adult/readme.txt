#*****Nov 25th, 2019*******

# run cmnd.txt 1-11 manually


#*****Jan 10th, 2020*******

# run cmnd.txt 12-14
qsub -q rg-el7 -N intersection.introns.exons.intergenic.adult -m bea -M beatrice.borsari@crg.eu -cwd -o /no_backup/rg/bborsari/projects/enhancers_neural_development/logs/$JOB_NAME -e /no_backup/rg/bborsari/projects/enhancers_neural_development/errors/$JOB_NAME cmnd.txt
