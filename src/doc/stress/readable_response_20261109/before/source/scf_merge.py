#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import os

def scf(case, wopt, outfile='.scf'):
    def SearchScf(ist,KEY,dat):
        for i in range(ist,len(dat)):
            line=dat[i].strip()
            if line[:len(KEY)]==KEY:
                break
        return (line,i)
    
    def GetSum(KEY, dat_up, dat_dn):
        line_up, iu_ = SearchScf(0,KEY,dat_up)
        line_dn, id_ = SearchScf(0,KEY,dat_dn)
        if iu_<len(dat_up) and id_<len(dat_dn):
            del dat_up[iu_]
            del dat_dn[id_]
            aSUM = 0.5*(float(line_up[30:50])+float(line_dn[30:50]))
            return line_up[:30]+("%20.9f"%aSUM)+line_up[50:]
        else:
            return ' '

    def GetAverageTensor(KEY, dat_up, dat_dn):
        """Remove two spin-resolved 3x3 tensors and return their average.

        In the antiferromagnetic DMFT workflow, dmft2 and dmft2dn each
        describe a complete, 62-electron MnTe cell with the two opposite
        assignments of the local self-energy.  Therefore every extensive
        valence quantity is averaged, rather than added.  GetSum applies this
        rule to the scalar Klein-functional terms and the force merger applies
        it to atomic forces.  Analytic stress introduced the three
        :STR_VAL001--003 rows, which need exactly the same treatment.

        The rows are removed from dat_up/dat_dn so that a downstream reader
        cannot accidentally select one spin mapping merely because it occurs
        last in the merged SCF stream.
        """
        rows_up = []
        rows_dn = []
        indices_up = []
        indices_dn = []

        for alpha in range(1, 4):
            tag = KEY + ('%03d:' % alpha)
            line_up, iu_ = SearchScf(0, tag, dat_up)
            line_dn, id_ = SearchScf(0, tag, dat_dn)

            # Stress is optional.  Return no rows unless both dmft2 channels
            # contain a complete tensor; silently mixing a partial tensor
            # would be much more dangerous than omitting it.
            if not line_up.startswith(tag) or not line_dn.startswith(tag):
                return []

            fields_up = line_up.split()
            fields_dn = line_dn.split()
            if len(fields_up) < 4 or len(fields_dn) < 4:
                return []

            rows_up.append([float(x) for x in fields_up[1:4]])
            rows_dn.append([float(x) for x in fields_dn[1:4]])
            indices_up.append(iu_)
            indices_dn.append(id_)

        # Delete in reverse order so the recorded indices remain valid.
        for index in sorted(indices_up, reverse=True):
            del dat_up[index]
        for index in sorted(indices_dn, reverse=True):
            del dat_dn[index]

        averaged = []
        for alpha in range(3):
            tag = KEY + ('%03d:' % (alpha + 1))
            values = [0.5 * (rows_up[alpha][beta] + rows_dn[alpha][beta])
                      for beta in range(3)]
            averaged.append(tag + '  ' + ''.join('%20.12f' % value
                                                   for value in values) + '\n')
        return averaged
    
    # First merger case.scf0, case.scf1 and case.scfso
    dat=[]
    scfs = ['0','1','so']
    if wopt['updn']: scfs = ['0','1up','1dn','so']
    for i in scfs:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()
    
    # case.scf2 is complicated when both up & dn are present
    # We need to average over case.scf2 and case.scf2dn
    if wopt['m_extn'] and os.path.exists(case+'.scf2'+wopt['m_extn']):
        lineo, io = SearchScf(0,':NATO',dat)
        nat = int(lineo.split(':')[2].split()[0])
        
        up__ = 'up' if wopt['updn'] else ''
        dat_up = open(case+'.scf2'+up__,'r').readlines()
        dat_dn = open(case+'.scf2'+wopt['m_extn'],'r').readlines()
        fu=case+'.scf2'+up__
        fd=case+'.scf2'+wopt['m_extn']

        # dmft2 and dmft2dn are the two complete antiferromagnetic spin
        # mappings.  Keep one physical valence-stress tensor, their mean.
        # GetAverageTensor also removes the two unaveraged copies from the
        # channel-specific sections appended below.
        lines_stress = GetAverageTensor(':STR_VAL', dat_up, dat_dn)
    
        line_s1 = GetSum(':SUM ', dat_up, dat_dn)
        line_s2 = GetSum(':XSUM ', dat_up, dat_dn)
        line_s3 = GetSum(':YSUM ', dat_up, dat_dn)
        line_s4 = GetSum(':ZSUM ', dat_up, dat_dn)
        line_s5 = GetSum(':WSUM ', dat_up, dat_dn)
        line_s6 = GetSum(':QSUM ', dat_up, dat_dn)
        line_s7 = GetSum(':VSUM ', dat_up, dat_dn)
        
        dat += '\n------- '+fu+' and '+fd+' --------\n'
        if lines_stress:
            dat += '\n------- averaged dmft2/dmft2dn valence stress -------\n'
            dat += lines_stress
    
        lineu, iu_ = SearchScf(0,':NOE ',dat_up)
        lined, id_ = SearchScf(0,':NOE ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[:iu_]
        dat += lineu
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[:id_]
        dat += lined
    
        lineu, _iu = SearchScf(iu_,':FER ',dat_up)
        lined, _id = SearchScf(id_,':FER ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += lineu
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        dat += lined
        iu_,id_ = _iu,_id
        
        lineu, _iu = SearchScf(iu_,':CHA ',dat_up)
        lined, _id = SearchScf(id_,':CHA ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        iu_,id_ = _iu,_id
        chup = float(lineu.split()[-1])
        chdn = float(lined.split()[-1])
        chtot = 0.5*(chup+chdn)
        dat += lineu[:40]+("%14.6f"%chtot)
    
    
        lineu, _iu = SearchScf(iu_,':DRHO',dat_up)
        lined, _id = SearchScf(id_,':DRHO',dat_dn)
            
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        iu_,id_ = _iu,_id
        
        drho_up = float(lineu.split()[1])
        drho_dn = float(lined.split()[1])
        #drho_up = float(lineu[23:39])        
        #drat_up = float(lineu[56:])
        #drho_dn = float(lined[23:39])        
        #drat_dn = float(lined[56:])
        drho = 0.5*(drho_up+drho_dn)
        #drat = 0.5*(drat_up+drat_dn)
        dat += lineu[:5]+("%20.12f "%drho)+'\n\n'
        if (abs(drho) > 0.5):
            print("WARNING : Electron charge density is very different than expected. You should change recomputeEF and set it to 1 (switch in on)")
            
        dat += line_s1+'\n'+line_s2+'\n'+line_s3+'\n'+line_s4+'\n'+line_s5+'\n'+line_s6+'\n'+line_s7+'\n\n'
        
    
        # merging forces
        lineu, iu_ = SearchScf(0,'VALENCE-FORCE',dat_up)
        lined, id_ = SearchScf(0,'VALENCE-FORCE',dat_dn)
        if iu_<len(dat_up)-1 and id_<len(dat_dn)-1:  # Forces are present
            for KEY in [':FVA',':FSU']:
                for iat in range(nat):
                    lineu, _iu = SearchScf(iu_,KEY,dat_up)
                    lined, _id = SearchScf(id_,KEY,dat_dn)
                    for i in range(iu_,_iu): # Printing what is between two forces :FVA/:FSU
                        dat += dat_up[i]
                    iu_ = _iu+1
                    id_ = _id+1
                    dlu = lineu.split()
                    dln = lined.split()
                    # average valence force between up and dn
                    Fa = [0.5*(float(dlu[i])+float(dln[i])) for i in range(3,6)]
                    # print out average valence force
                    dat += dlu[0]+('%9s'%dlu[1])+('%17.9f'*4) % tuple([sqrt(dot(Fa,Fa))]+Fa),'\n'
    else:
        # for non-magnetic calculation it is very simple
        if os.path.exists(case+'.scf2'):
            dat += '\n------- '+case+'.scf2'+' --------\n'
            dat += open(case+'.scf2','r').readlines()
            
    if os.path.exists('Eorb_imp.dat'):
        dat += '\n---------- Eorb_imp.dat --------------\n'
        dat += open('Eorb_imp.dat','r').readlines()
        
    scfs = ['1s','2s','c']
    if wopt['updn']: scfs = ['cup','cdn']
    for i in scfs:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()

    # Publish only after valence, density and fresh core parts are assembled.
    # The helper replaces their total and verifies the exact response snapshot.
    from stress_response_publish import response_rows, mark_published
    correction = response_rows(case, ''.join(dat))
    dat += correction
    with open(case+outfile, 'a') as fs:
        fs.writelines(dat)
    mark_published(case, correction, case+outfile)

def scfm(fday, case, WIEN):
    dat=[]
    for i in ['m']:
        dat += '\n------- '+case+'.scf'+i+' --------\n'
        if os.path.exists(case+'.scf'+i):
            dat += open(case+'.scf'+i,'r').readlines()
            
    with open(case+'.scf', 'a') as fs:
        fs.writelines(dat)

if __name__ == '__main__':
    import os
    import utils
    w2k = utils.W2kEnvironment()    # W2k filenames and paths
    inldn = None
    m_extn = 'dn' if os.path.exists(w2k.case+'.indmfl'+'dn') else ''
    wopt={'m_extn':m_extn, 'updn':True}
    scf(w2k.case, wopt, '.scf_test')
    
