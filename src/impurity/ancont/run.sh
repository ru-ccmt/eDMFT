#!/bin/bash


grep -v '#' Sig.average > Sig.work

strans Sig.work '#1' '#2'  '#3'   '#12' '#13' > Sig.1
strans Sig.work '#1' '#4'  '#5'   '#14' '#15' > Sig.2
strans Sig.work '#1' '#6'  '#7'   '#16' '#17' > Sig.3
strans Sig.work '#1' '#8'  '#9'   '#18' '#19' > Sig.4
strans Sig.work '#1' '#10' '#11'  '#20' '#21' > Sig.5

for x in `seq 1 5`; do mkdir n.$x; done

for x in `seq 1 5`; do cp Sig.$x pmesh.py ancont.py chi2f.so n.$x; done



for x in `seq 1 5`; do cd n.$x; ./ancont.py -sig Sig.$x -nom 100 -beta 100. -wexp 1.15 -Ng 60 -FermiLiquid True -L0 20 -lcut 0.18 -b 0.9 -alpha4 1. -p0 0.1 -vunity 0.01 ; cd ../; done



strans n.1/Sig.out n.2/Sig.out n.3/Sig.out n.4/Sig.out n.5/Sig.out '#1' '#1:2' '#1:3' '#2:2' '#2:3' '#3:2' '#3:3' '#4:2' '#4:3' '#5:2' '#5:3'  > sig.inp

#cp sig.inp onreal

#cd onreal; ./iterate.py ham.dat; cd ../
