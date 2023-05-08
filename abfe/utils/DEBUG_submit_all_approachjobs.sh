#!/bin/sh

initPATH=${PWD}

for x in $(ls -d ${initPATH}/*/*); 
do 
  cd ${x}
  schedulerPath="${x}/scheduler.sh"
  echo "Submit: ${schedulerPath}";
  snakemake --unlock;
  cd complex
  snakemake --unlock;
  cd ..
  cd ligand
  snakemake --unlock;
  cd ..
  eval ${schedulerPath}; 
  cd ${initPATH}
done
