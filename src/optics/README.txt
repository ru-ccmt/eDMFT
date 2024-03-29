This is version of the code from 2014, which has an option to use
finite temperature Fermi function broadening.
-----------------------------------------------
1) Run DMFT
2) Obtain matrix elements:

   Prepare "case.inop" file. Check Wien2K manual.
   
   in case.inop change MME to ON
   > x optic[c] [-so] [-up]

   You should get "case.mommat" and "case.symop"
   
   In the new version of wien2k you might need to rename case.mommat2 to case.mommat, 
   or just create a link
   ln -s case.mommat2 case.mommat

3) run DMFT in mode "u":
   a) execute
   > x_dmft.py dmftu
   You should get two files: Udmft[updn].0 and BasicArrays.dat

4) prepare input file named "dmftopt.in", which contains:

    0          # Temperature  (this is new in 2014)
    sr2ruo4    # case
    0.01       # gamma  -- broadening of all bands / we recommend to keep nonzero value
    0.0        # gammac -- broadening of the correlated bands (in addition to the self-energy)
    4          # ommax  -- maximum frequency for the optics calculation
    1e-2       # delta  -- minimum separation of frequency for logarithmic mesh in frequency
    5          # Nd     -- number of points in each linear mesh, a subset of logarithmic mesh
    F          # Qsym: [F|T]. Do we need to symmetrize? (over all k-points or just irreducible.) If Qsym is F, it goes over irreducible, if Qsym is T, it goes over all.
    F          # InterbandOnly [F|T] (F -- all, T--interband)
    10         # dwindow -- Not all bands are used in computation of optics, but only those from [-omeg-dwindow, omega+dwindow]. You should increase this to make sure you are not cutting some bands which are needed.
    2          # Ndirection -- How many optics type do you want to compute: xx, yy, zz,...
    0.5 0.0 0.0    # alphaV(1,:)
    0.0 0.5 0.0    # alphaV(2,:)
    0.0 0.0 0.0    # alphaV(3,:)
    0.0 0.0 0.0    # alphaV(1,:)
    0.0 0.0 0.0    # alphaV(2,:)
    0.0 0.0 1.0    # alphaV(3,:)


   prepare self-energy sig.inp1 by
   > ssplit.py

   and finally run DMFT optics by
   > dmftopt
   
   You will need "case.energy" , "case.mommat", "case.symop" ,
   "Udmft.0" and "BasicArrays.dat".

   Result is stored in "optics.dat". The corresponding total density
   of states is in "optdos.dat".

----------------------------------------------
To perform finite temperature calculation, you just set T in the first line of dmftopt.in file.
To speed up the calculation, you should choose ommax to be around 5-times the value of temperature. 
Too large value of ommax will result in very slow calculation, because we use 
linear mesh with the first point at delta, and last point at ommax. The input value of Nd is ignored,
and the number of points is choosen such that the mesh  x=[-ommax,......,-delta,0,delta,.....ommax] is linear.
The mesh x is used in the integration. The frequency mesh for sigma(w) will have even less points, because
when we use cutoof ommax in integration over x, sigma(w) can not be computed on all points.

In the output, you will see the the number of points 'Nw=...', which is the number of points in the x-linear mesh.
We do not use second or third level of logarithmic mesh, hence N0=1.
You can also check on how many points sigma(w) will be calculated. Check 'istart=1 iend=...'. The value of iend*delta
is the largest frequecy calculated.
