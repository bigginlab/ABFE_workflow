#!/usr/bin/env bash

ori=${PWD}
snakemake --unlock -c1
snakemake -R $(snakemake --list-code-changes) --touch -c4

for dir in */?
do
  cd ${dir}
  snakemake --unlock -c1
  snakemake -R $(snakemake --list-code-changes) --touch -c4
  cd ${ori}
done

for dir in */?/ligand
do
  cd ${dir}
  snakemake --unlock -c1
  snakemake -R $(snakemake --list-code-changes) --touch -c4
  cd ${ori}
done

for dir in */?/complex
do
  cd ${dir}
  snakemake --unlock -c1
  snakemake -R $(snakemake --list-code-changes) --touch -c4
  cd ${ori}
done
