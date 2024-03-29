#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import optparse
#from scipy import *
import sys, re
import os
import glob

def TOTtoFOR(case):
    case_in2s = glob.glob(case+'.in2*')
    for fin2 in case_in2s:
        mode = open(fin2,'r').readline().split()[0]
        if mode!='FOR':
            fw = open(fin2,'r')
            data = fw.readlines()   # reads the entire file
            fw.close()              # needs to close to safely reopen
            data[0] = 'FOR'+data[0][3:] # changing the first line
            fw = open(fin2,'w')     # now rewriting it
            fw.writelines(data)     # writes data back to file
            fw.close()              # safely closing it.
        
def AreWeRelaxingStructure(case):
    mix_method = open(case+'.inm','r').readline().split()[0]
    return (mix_method=='MSR1a' or mix_method=='MSR2a' or mix_method=='MSECa')

def AreWeRelaxingStructure2(case):
    mix_method = open(case+'.inm','r').readline().split()[0]
    RelaxingStructure = (mix_method=='MSR1a' or mix_method=='MSR2a' or mix_method=='MSECa')
    return (RelaxingStructure, mix_method)

def RelaxAgain(case, mix_method):
    case_inm = case+'.inm'
    fw = open(case_inm,'r')
    data = fw.readlines()   # reads the entire file 
    fw.close()              # needs to close to safely reopen 
    old_mix_method = data[0].split()[0]
    data[0] = mix_method + data[0][len(mix_method):]
    fw = open(case_inm,'w') # now rewriting it    
    fw.writelines(data)     # writes data back to file 
    fw.close()              # safely closing it.

def AreForcesConverged(case, how_many_times=3, min_Ene_diff=0.001, min_Rho_diff=0.01, info=sys.stdout):
    scf_content = open(case+'.scf').readlines()

    line_index=[]
    for i0 in range(len(scf_content)-1,0,-1):
        if scf_content[i0][:5]==':FRMS':  # Finds the force condition
            line_index.append(i0)         # and remembers the line number
            if len(line_index)>=how_many_times: break  # We need only 3 occurences
    
    Differences=[]                        # This will contain pairs of energy and charge differences for the last three iterations
    for ii in line_index:                 # ii is the line with force condition
        lined = scf_content[ii].split()
        if len(lined)>=13 :               # Forces are being computed if there are 13 items
            ene=[]                        # Now checking how well is energy and charge converged
            rho=[]
            for j in range(ii-1,0,-1):    # Reading the scf file backward from the force condition
                if scf_content[j][:4]==':ENE':    # save any energy we find backward
                    ene.append( float(scf_content[j].split()[-1]) ) # energy at this and previous iteration
                if scf_content[j][:4]==':DIS':    # save any charge we find backaward
                    rho.append( float(scf_content[j].split()[-1]) ) # charge at this and previous iteration
                if len(ene)>=2 and len(rho)>=2:          # Need only two energies and charges
                    break
            e_diff=1.
            if len(ene)>=2: e_diff = abs(ene[1]-ene[0])
            r_diff=1.
            if len(rho)>=1: r_diff = abs(rho[0])
            Differences.append((lined[12],e_diff, r_diff))
    
    for item in Differences:
        print('Convergence of Force=', item[0], ', energy=', item[1], '('+str(item[1]<min_Ene_diff)+')',' and charge=', item[2], '('+str(item[2]<min_Rho_diff)+')', file=info)

    success=False
    if len(Differences)>=how_many_times:
        success=True
        for i in range(how_many_times):
            success = success and (Differences[i][0]=='T' and Differences[i][1]<=min_Ene_diff and Differences[i][2]<=min_Rho_diff)
    if success:
        print('Structure is now converged', file=info)
        
    return success


def StopStructureOptimization(case,info=sys.stdout):
    fw = open(case+'.inm','r')
    data = fw.readlines()   # reads the entire case.inm file
    fw.close()              # needs to close to safely reopen
    data[0] = re.sub('MSR1a', 'MSR1 ', data[0])    # change MSR1a to MSR1
    data[0] = re.sub('MSECa', 'MSEC3', data[0])    # change MSECa to MSEC3
    data[3] = re.sub('(9999\s*)8', '\\1 10', data[3]) # change 8 to 10 in history
    fw = open(case+'.inm','w')     # now rewriting it
    fw.writelines(data)     # writes data back to file
    fw.close()              # safely closing it.
    
    print('MSR1a/MSECa changed to MSR1/MSEC3 in '+case+'.inm, relaxing only the electronic structure.', file=info)
    
    # remove broydens
    for filename in glob.glob("*.broy*"):
        os.remove(filename)
    print('....broydens removed', file=info)
    
    dscfm = open(case+'.scfm','r').readlines() # now reading the entire case.scfm file
    dscfm = dscfm[::-1]  # all lines are swapped, so that we can now search from the beginning
    greed=None
    for line in dscfm:
        m = re.search('GREED', line)
        if m is not None:
            dat = line.split()
            for i in range(len(dat)):
                if  dat[i][:5]=='GREED': 
                    greed = float(dat[i+1])
                    break
        if greed is not None:
            break

    if greed is not None:
        new_greed = max(greed/2,0.05)
        fw = open('.msec', 'w')
        print('%5.3f' % (new_greed,), file=fw)
        fw.close()
        print('At the end of structural optimization the GREED was %5.3f and is now set to %5.3f' % (greed, new_greed), file=info)
    
if __name__ == '__main__':
    import utils
    env = utils.W2kEnvironment()
    case = env.case

    (RelaxingStructure, mix_method) = AreWeRelaxingStructure2(case)
    print('mix_method = ', mix_method)
    #RelaxAgain(case, mix_method+'a')
    
    #succ = AreForcesConverged(case, how_many_times=3)
    #print('Structure converged =', succ)
    
    sys.exit(0)
    
    min_Ene_diff = 0.001
    min_Rho_diff = 0.01
    how_many_times=3

    testmsr = AreWeRelaxingStructure(case)
    if testmsr:
        TOTtoFOR(case)
    

    succ = AreForcesConverged(case, how_many_times=3)
    print('Structure converged =', succ)
    
    StopStructureOptimization(case)
