#!/bin/bash

echo 0 | gmx trjconv -s ../npt_prod1/npt_prod1.tpr -f ../npt_prod1/npt_prod1.xtc -o npt_prod1_whole.xtc -pbc whole

echo 0 | gmx trjconv -s ../npt_prod1/npt_prod1.tpr -f npt_prod1_whole.xtc -o npt_prod1_NOJUMP.xtc -pbc nojump

gmx trjconv -s ../npt_prod1/npt_prod1.tpr -f npt_prod1_NOJUMP.xtc -o npt_prod1_center -pbc mol -center -ur compact << EOF
1
0
EOF

rm npt_prod1_whole.xtc npt_prod1_NOJUMP.xtc
