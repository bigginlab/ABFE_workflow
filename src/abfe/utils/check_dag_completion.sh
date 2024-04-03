#!/usr/bin/env bash

ori=${PWD}
snakemake --unlock -c1 -s Snakefile.smk
snakemake -R $(snakemake --list-code-changes -s Snakefile.smk) --touch -c4 -s Snakefile.smk

for dir in */?
do
  cd ${dir}
  snakemake --unlock -c1 -s Snakefile.smk
  snakemake -R $(snakemake --list-code-changes -s Snakefile.smk) --touch -c4 -s Snakefile.smk
  cd ${ori}
done

for dir in */?/ligand
do
  cd ${dir}
  snakemake --unlock -c1 -s Snakefile.smk
  snakemake -R $(snakemake --list-code-changes -s Snakefile.smk) --touch -c4 -s Snakefile.smk
  cd ${ori}
done

for dir in */?/complex
do
  cd ${dir}
  snakemake --unlock -c1 -s Snakefile.smk
  snakemake -R $(snakemake --list-code-changes -s Snakefile.smk) --touch -c4 -s Snakefile.smk
  cd ${ori}
done
