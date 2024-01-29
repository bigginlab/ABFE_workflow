#!/bin/sh

initPATH=${PWD}

for x in $(ls -d ${initPATH}/*/*); 
do 
  cd ${x}
  schedulerPath="${x}/scheduler.sh"
  echo "Submit: ${schedulerPath}";
  snakemake --unlock -s Snakefile.smk;
  cd complex
  snakemake --unlock -s Snakefile.smk;
  cd ..
  cd ligand
  snakemake --unlock -s Snakefile.smk;
  cd ..
  #eval ${schedulerPath}; 
  cd ${initPATH}
done
