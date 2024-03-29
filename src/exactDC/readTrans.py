from numpy import array, zeros
# @Copyright 2007 Kristjan Haule

def ReadTrans(filename, fh_info):
    """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
    fh = open(filename, 'r')
    data = fh.readlines()
    
    (n1,n2) = list(map(int, data[0].split()[:2]))

    Sigind=[]
    for i in range(n1):
        Sigind.append( list(map(int, data[i+2].split()[:n2])) )
    Sigind = array(Sigind)

    #print >> fh_info, 'len(data)', len(data)
    #print >> fh_info, 'n1=', n1
    if len(data) >= n1+n1+3:
        n2 = n1
        CF=[]
        for i in range(n2):
            cl = array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
    elif len(data)>=n1+n1/2+3:
        n2 = n1/2
        CF=[]
        for i in range(n2):
            cl = array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
        CFN = zeros((2*n2,2*n2), dtype=complex)
        CFN[:n2,:n2] = CF
        CFN[n2:,n2:] = CF
        CF = CFN
    else:
        CF = identify(n1)
    #print >> fh_info, 'CF=', CF
    
    return (Sigind, CF)

if __name__ == '__main__':
    import sys
    (Sigind, CF) = ReadTrans('Trans.dat', sys.stdout)
    print('Sigind=', Sigind)
    print('CF=', CF)
