import os
from random import seed
from random import randint
for i in range(0,2):
    f = open('command.sh','w')
    g4seed = randint(0,10000000)
    f.write('#!/bin/bash -l \n')
    f.write('#$ -l h_rt=12:00:0 \n')
    f.write('#$ -l mem=5G \n')
    f.write('#S -wd /home/rmapcco/Software/ProtonMC \n')
    f.write('source /home/rmapcco/Software/root/bin/thisroot.sh \n')
    f.write('source /home/rmapcco/Software/geant4_10_01_p03/share/Geant4-10.1.3/geant4make/geant4make.sh \n')
    f.write('cd /home/rmapcco/Software/ScintillatorMC \n')
    f.write('/home/rmapcco/Software/ScintillatorMC/bin/Linux-icc/ScintillatorMC 1000000 200 LasVegas 0 30 ' + str(g4seed) + ' 1')
    f.close()
    os.system('qsub command.sh')
    os.system('rm command.sh')
