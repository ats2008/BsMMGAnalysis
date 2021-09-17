#!/usr/bin/env python3
import sys
import os

if len(sys.argv)>1:

    fname=sys.argv[1]
    print('coustomizing ',fname)
    f=open(fname)
    fout=open(fname+'.new','w')
    l=f.readline()
    while l:
        if '_customInfo' in l:
            l='#'+l
        fout.write(l)
        l=f.readline()
    
    fout.write('\n')
    fout.write('# Bs2MMGCustomization\n')
    fout.write('from customizeForBs2MMG import *\n')
    fout.write('process = customizeBsToMMGTrigPathForEfficiencyFiles(process)\n\n')
    #fout.write('from HLTrigger.Configuration.customizeHLTforCMSSW import customiseFor2018Input \n')
    #fout.write('customiseFor2018Input(process) \n\n')
    fout.write('#Customization for EG Debug\n')
    fout.write('from HLTrigger.Configuration.customizeHLTforEGamma import * \n')
    fout.write('process = customiseEGammaMenuDev(process)\n\n')
    fout.write('process.options.numberOfThreads=1')
    
    f.close()
    fout.close()
    os.system('mv '+fname+' '+fname+'.bak')
    os.system('mv '+fname+'.new'+' '+fname)
else:
    print("usage ./putcust...py CFG_FILENAME")
