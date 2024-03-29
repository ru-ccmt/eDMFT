1) run Wien2K and prepare case.indmfl
2) run kgena.py to generate all k-points in 1st-BZ
  >kgena.py "[5,5,5]"
3) run lapw1 on new klist
  >x lapw1
4) print projector Udmft.0
 >x_dmft.py dmftu -g
  # Here  -g stands for turning off the loop over group operations. 
  # We use all k-points and hence we do not go over all group operations
5) execute downfold.py
 >downfold.py -E -2.0 -n 12 
  # Here -E stands for the lower cuttof for the bands (in eV from EF) and -n 
  # contains number of bands considered in downfolding.

You will obtain "hamiltonian.dat"

Use iterate.py to plot band or perform DMFT with the tight binding hamiltonian


