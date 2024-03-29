import re,sys,os
# @Copyright 2007 Kristjan Haule
def FindChemicalPotential(case, updn):
    # Looking for the LDA chemical potential
    Ry2eV = 13.60569193
    EF_found = False

    if os.path.isfile(case+'.in2'):
        fname = case+".in2"
    elif os.path.isfile(case+'.in2c'):
        fname = case+".in2c"
    else:
        raise Exception("Failed to determine number of electrons. You should have case.in2 or case.in2c file!")
        
    lines = open(fname, 'r').readlines()
    NOE = float(lines[1].split()[1])
    print('NOE=', NOE)
        
    fname = case+".scf2"
    if os.path.isfile(fname):
        fscf = open(fname, 'r')
        print("opening", fname)
    elif os.path.isfile(fname+updn):
        fscf = open(fname+updn, 'r')
        print("opening", fname+updn)
    else:
        fscf = None

    if fscf is not None:
        lines = fscf.readlines()
        for line in lines:
            if re.match(r':FER', line) is not None:
                EF = float(line[38:])*Ry2eV
                print('Ef=', EF)
                EF_found = True

    # The previous DMFT chemical potential
    if  os.path.isfile('EF.dat'):
        fmu = open('EF.dat','r')
        mu = float(next(fmu))
        print('Found DMFT-mu=', mu)
        EF = mu
        EF_found = True

    if not EF_found:
        fname = case+".scf"
        if os.path.isfile(fname):
            print("opening", fname)
            f = open(fname, 'r')
            lines = f.readlines()
            if not EF_found:
                for line in lines:
                    if re.match(r':FER', line) is not None:
                        EF = float(line[38:])*Ry2eV
                        print('Ef=', EF)
                        EF_found = True

    if not EF_found:
        raise Exception("Failed to determine chemical potential.")

    return (EF, NOE)

if __name__ == '__main__':
    import utils
    w2k = utils.W2kEnvironment()
    updn=''
    EF,NOE = FindChemicalPotential(w2k.case, updn)
    print('EF=', EF, 'NOE=', NOE)
