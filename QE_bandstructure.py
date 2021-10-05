#!/usr/bin/env python2

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib import gridspec
import argparse
from argparse import RawTextHelpFormatter
from itertools import islice
import re
import os

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', action='store', dest='inpt',
                    help='output of bands calculation', required=True, type=str)
parser.add_argument('-fermi', action='store', dest='fermi',
                    help='put "ins" for insulator and its value for metals', required=True)
parser.add_argument('-emin', action='store', dest='emin',
                    help='emin', required=False, type=float, default=-4)
parser.add_argument('-emax', action='store', dest='emax',
                    help='emax', required=False, type=float, default=4)

results = parser.parse_args()

###########################################################
# Read the name of the high-symmetry points
bandfile = results.inpt
namefile = bandfile.strip().split(".")[0]
Emin = results.emin
Emax = results.emax

kpointfile = open('./KPOINTS', 'r')
kpointfile.readline()
nump_of_points = int(kpointfile.readline().split()[0])
kpointfile.readline()
kpointfile.readline()
contents = kpointfile.read().split('\n')

kname = []
for i in range(0, len(contents), 3):
    kname.append(contents[i].split()[-1].strip('!''#'))

kpointfile.seek(0)
for line in kpointfile:
    pass
last = line

lsym = last.split()[-1].strip('!''#')  # Last_Kpoint
kname.append(lsym)  # Just we need kname from this part

for i, item in enumerate(kname):
    if item == 'G':
        kname[i] = '$\Gamma$'
    elif item == '\Gamma':
        kname[i] = '$\Gamma$'
        # print klist
###########################################################
flist = open(bandfile).readlines()

numelec = []
parsing = False
cnt = 0
eigen_band = []
kp = 0
for line in flist:
    if "number of electrons" in line:
        el = float(line.strip().split()[4])
        numelec.append(el)
    if "  k =" in line:
        cnt = cnt + 1
        parsing = True
    elif "  k =" in line:
        parsing = False
    if parsing:
        if cnt <= 1:
            if not "  k =" in line:
                linefloat = [float(i) for i in line.strip(
                ).strip().replace('-', ' -').strip().split()]
                if linefloat != []:
                    #                    print linefloat
                    eigen_band.append(linefloat)

numelec = numelec[0]
nk = int(cnt)
nbnd = int(numelec/2)
eigen_band2 = np.array(eigen_band)
lnum = len(eigen_band2)
nband = 0
for n in range(eigen_band2.shape[0]):
    nband += len(eigen_band2[n])
eigen_all = []

kvector = []
j = 0
record = False
count = 0
with open(bandfile) as f:
    for line in f:
        if "  k =" in line:
            record = True
        elif "Writing output data" in line:
            count = count+1
            record = False
        if record:
            if count <= 1:
                if "  k =" in line:
                    kvec = line.replace('  k =', '').replace(
                        '-', ' -').strip().split()[0:3]
                    kvec00 = map(float, kvec)
                    kvector.append(kvec00)
                    j = j + 1
                    b = np.empty((0))
                    for i in range(lnum+1):
                        a = f.next().strip().replace('-', ' -').split()
                        a = np.array(a, dtype=float)
                        eigen_all.extend(a)
                        # if a!=[]:                        #This part makes problem
                        #     for t in range(len(a)):      #this part makes problem
                        #         eigen_all.append(float(a[t]))
eigen_all2 = np.array(eigen_all)
eigen_all3 = np.reshape(eigen_all2, nk*nband)
kvector = np.array(kvector)


dist = 0
kdist = []
kdist.append(0.0)
for i in range(len(kvector)):
    for j in range(len(kvector)):
        if j-i == 1:
            dist += np.linalg.norm(kvector[j]-kvector[i])
            kdist.append(dist)
kdist = np.array(kdist)

kdist_hsym = []
for i in range(0, len(kdist), nump_of_points-1):
    kdist_hsym.append(kdist[i])


Eigentotal = []
for j in range(nband):
    eigen_per_kpt = []
    for i in range(j, len(eigen_all3), nband):
        eigenext = eigen_all3[i]
        eigen_per_kpt.append(eigenext)
    Eigentotal.append(eigen_per_kpt)

Eigentotal = np.array(Eigentotal)


if results.fermi == "ins":
    Ef = float(np.amax(Eigentotal[(int(numelec/2)-1)]))
    BandGap = abs(Ef-float(np.amin(Eigentotal[(int(numelec/2))])))
    HOMO = []
    for i in range(nbnd):
        for j in range(nk):
            if Eigentotal[i][j] <= Ef:
                HOMO.append(i)
    HOMO = np.amax(HOMO)
    TVB = np.amax(Eigentotal[HOMO])
    BCB = np.amin(Eigentotal[HOMO+1])
    gap = abs(TVB-BCB)
    middle_point = np.amax(Eigentotal[HOMO])+gap/2.0
    Fermi_dist_to_middle = Ef-middle_point
else:
    Ef = float(results.fermi)
    HOMO = []
    for i in range(nbnd):
        for j in range(nk):
            if Eigentotal[i][j] <= Ef:
                HOMO.append(i)
    HOMO = np.amax(HOMO)
    TVB = np.amax(Eigentotal[HOMO])
    BCB = np.amin(Eigentotal[HOMO+1])
    gap = abs(TVB-BCB)
    middle_point = np.amax(Eigentotal[HOMO])+gap/2.0
    Fermi_dist_to_middle = Ef-middle_point

print "----Check if VB and CB are the same with those of SCF-----"
print "Top of the VB       = ", TVB, "eV"
print "Bottom of  CB       = ", BCB, "eV"
print "Band  Gap  is       =  ", round(gap, 4), "eV"
print "Fermi Energy Ef     = ", Ef
print "Band Gap Center BGC = ", middle_point
print "HOMO Index          = ", HOMO+1
print "LUMO Index          = ", HOMO+2
print "----------------------------------------------------------"


font_manager.findfont('Helvetica world')
plt.rc('font', family='Helvetica world')
plt.rc('font', serif='Helvetica world', size=11)
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['xtick.minor.width'] = 1.0
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['ytick.minor.width'] = 1.0
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'

fig = plt.figure(figsize=(5, 3))
gs = gridspec.GridSpec(1, 1)
ax0 = fig.add_subplot(111)
ax1 = fig.add_subplot(gs[0])
ax0.spines['top'].set_color('none')
ax0.spines['bottom'].set_color('none')
ax0.spines['left'].set_color('none')
ax0.spines['right'].set_color('none')
ax0.set_xticks([])
ax0.set_yticks([])

#
for i in range(nband):
    m = Eigentotal[i]-Ef
    ax1.plot(kdist, m, lw=1, color='navy')

fout = open('BANDGNUEfzero.dat', 'w')
for j in range(nband):
    for k in range(len(kdist)):
        fout.write("%12d %18.8f %18.8f\n" %
                   (k+1, kdist[k], Eigentotal[j][k]))
    fout.write("\n")

fouthomo = open('BAND_HOMO.dat', 'w')
fouthomo_1 = open('BAND_HOMO_1.dat', 'w')
foutlumo = open('BAND_LUMO.dat', 'w')
foutlumo_1 = open('BAND_LUMO_1.dat', 'w')
for j in range(nband):
    m = np.array(Eigentotal[j])
    if j == int(HOMO):
        for i in range(len(kdist)):
            fouthomo.write("%8d %18.8f %18.8f\n" % (i+1, kdist[i], m[i]))
    if j == int(HOMO-1):
        for i in range(len(kdist)):
            fouthomo_1.write("%8d %18.8f %18.8f\n" % (i+1, kdist[i], m[i]))
    if j == int(HOMO+1):
        for i in range(len(kdist)):
            foutlumo.write("%8d %18.8f %18.8f\n" % (i+1, kdist[i], m[i]))
    if j == int(HOMO+2):
        for i in range(len(kdist)):
            foutlumo_1.write("%8d %18.8f %18.8f\n" % (i+1, kdist[i], m[i]))


ax1.axis([0, np.amax(kdist), -5, 5])
ax1.axhline(0, color='red', lw=1, ls='-.')
ax1.set_ylabel("$E-E_{\mathrm{F}}$ [eV]")
ax1.set_xlabel(r"$\vec{k}$")

cwd = os.getcwd()
nameon = os.path.basename(cwd).split('-')[0]
bbox_props = dict(boxstyle="round", fc="khaki", ec="0.5", alpha=1, pad=0.2)
ax1.annotate(nameon, xy=(0.1, Emax-0.2*Emax), xycoords='data',
             horizontalalignment='left', verticalalignment='bottom', bbox=bbox_props, fontsize=12)


ax1.set_ylim(Emin, Emax)

for i in range(len(kdist_hsym)):
    ax1.axvline(kdist_hsym[i], -5., 5., color='black',
                lw=1, ls='--', alpha=0.5)
plt.xticks(kdist_hsym, kname, family='Serif')

plt.savefig(str(namefile)+'.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('Fig.png', bbox_inches='tight', pad_inches=0.02, dpi=1000)
# plt.show()