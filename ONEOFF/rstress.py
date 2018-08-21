#!/usr/local/bin/python
#python script to read in and plot data from A.C. Boley's rstress.f90 program.

import numpy
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys
import os
 
def main(*args):
    
    time1  = 12.0
    time2  = 18.0
    colors = ['k','b','g','r']
    istart = []
    iend   = []
    iskip  = []
    prefix = []
    pi = 3.14159265358979323846
    istart.append(205000)
    iend.append(500000)
    iskip.append(5000)
    prefix.append('/Users/scamicha/Research/Data/INITCOND/RES_STUDY/RSTRESS/rs_2cell_patches_LMAX64.')
    istart.append(110007)
    iend.append(385007)
    iskip.append(5000)
    prefix.append('/Users/scamicha/Research/Data/INITCOND/RES_STUDY/RSTRESS/rs_2cell_patches_LMAX128.')
    istart.append(205000)
    iend.append(650000)
    iskip.append(5000)
    prefix.append('/Users/scamicha/Research/Data/INITCOND/RES_STUDY/RSTRESS/rs_2cell_patches_LMAX256.')
    istart.append(205000)
    iend.append(980000)
    iskip.append(5000)
    prefix.append('/Users/scamicha/Research/Data/INITCOND/RES_STUDY/RSTRESS/rs_2cell_patches_LMAX512.')

    jmax    = 512
    for k in range(4):
        radius  = []
        rstress = []
        count   = 0
        i = istart[k]
        while i <= iend[k]:
            filename = prefix[k]+str(i)
            if not os.path.isfile(filename):
                i += iskip[k]
                continue
            f = open(filename,'r')
            teststr = f.readline()
            testarr = teststr.split()
            if float(testarr[0]) < time1 or float(testarr[0]) > time2:
                i += iskip[k]
                continue
            count += 1
            if count == 1:
                for line in f:
                    entries = line.split()
                    radius.append(float(entries[1]))
                    rstress.append(abs(float(entries[2])))
            else:
                for j,line in enumerate(f):
                    entries = line.split()
                    rstress[j] += abs(float(entries[2]))
            
            i += iskip[k]
        dz = radius[1]-radius[0]

        for i,value in enumerate(rstress):
            rstress[i] = value*2.0*pi*dz*pow(radius[i],2.0)/(count*64*(k+1))
        
        if k == 0:
            xmin = 10.0
            ymin = 0.0
            xmax = 50.0
            ymax = 7.0e38u
            axay = [xmin,xmax,ymin,ymax]
            print axay
            plt.axis(axay)
        plt.plot(radius,rstress,colors[k])
        
        del radius
        del rstress

    plt.savefig('2cell_patches_comp_abs.eps')
#    plt.show()
    
    return 0

if __name__ == "__main__":
    main()
