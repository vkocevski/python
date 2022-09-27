#! /usr/bin/env python
#
# Script that generates supercell within a user defined range of
# number of atoms and finds the supercell with smallest difference
# in lattice parameters a, b and c.
#
# It runs mcsqs to generate possible suppercels in sqscell.out.
# Then trims the sqscell.out to user defined number of cells having
# very similar lattice parameters and user defined anlges.
# It prints out the trimed cells in ascending order of lattice
# parameter error, with an error below user defined tolerance.

import numpy as np
import argparse
import os

# Options
parser = argparse.ArgumentParser()
parser.add_argument('--nat', type=int, required=True, help='Number of atoms in unit cell')
parser.add_argument('--nsqs', type=int, required=True, help='Number of supercells to be tested')
parser.add_argument('--nstart', type=int, required=True, help='Number of starting supercell')
parser.add_argument('--ang', type=float, required=True, help='Lowest angle of the supercell (in degrees')
parser.add_argument('--da', type=float, required=True, help='Lattice parameter tolerance (in Ang)')
parser.add_argument('--nc', type=int, required=True, help='Numer of cells')
args = parser.parse_args()

nat = args.nat
nstart = args.nstart
nsqs = args.nsqs+nstart
lima1 = args.ang
lima2 = 180-lima1
da = args.da
ncell = args.nc

# Run mcsqs to get lattice vectors for all possible supercells
def main():
    for ns in range(nstart-1,nsqs-1):
        os.system("timeout 2m ./mcsqs -n %s" % str((ns+1)*nat))
        trim_cell(ns)

# Triming sqscell.out
def trim_cell(ns):
    with open('sqscell.out') as f1:
        lines = f1.readlines()
    tot = int(lines[0].split()[0])  # total number of cells
    count = 0
    arr1 = np.zeros((ncell,3,3)) # initializing array to store trimed cell
    dlp_min = []
    dlp_min = [1000 for i in range(ncell)] # initializing array to store lattice parameter error
    for i in range(tot):
        a=[0, 0, 0]
        b=[0, 0, 0]
        c=[0, 0, 0]
        lp=[0, 0, 0]
        dt=[0, 0, 0]
        ang=[0, 0, 0]
        dlp = [0, 0, 0]
        # Reading sqscell.out
        a = [ float(x) for x in lines[4*i+2].split() ]
        b = [ float(x) for x in lines[4*i+3].split() ]
        c = [ float(x) for x in lines[4*i+4].split() ]
        # Calculating lattice parameters
        lp[0] = np.sqrt(np.sum(np.power(a,2)))
        lp[1] = np.sqrt(np.sum(np.power(b,2)))
        lp[2] = np.sqrt(np.sum(np.power(c,2)))
        # Calculating lattice angles
        dt[0] = np.sum(np.multiply(a,b))
        dt[1] = np.sum(np.multiply(a,c))
        dt[2] = np.sum(np.multiply(b,c))
        ang[2] = np.arccos(dt[0]/(lp[0]*lp[1]))*180/np.pi
        ang[1] = np.arccos(dt[1]/(lp[0]*lp[2]))*180/np.pi
        ang[0] = np.arccos(dt[2]/(lp[1]*lp[2]))*180/np.pi
        # Checking if lattice angles are withing user defined angles
        if ang[0] >= lima1 and ang[0] <= lima2 and ang[1] >= lima1 and ang[1] <= lima2 and ang[2] >= lima1 and ang[2] <= lima2:
            # Calculating difference in lattice parameters
            dlp[0] = abs(lp[0] - lp[1])
            dlp[1] = abs(lp[0] - lp[2])
            dlp[2] = abs(lp[1] - lp[2])
            # Checking if lattice parameter difference is below user defined tolerance
            if (dlp[0] < da and dlp[1] < da and dlp[2] < da):
                dlp2 = np.sum(np.power(dlp,2)) # error square
                count += 1 # counting how many cells satisfy the user defined parameters
                for n in range(ncell):
                    if dlp2 < dlp_min[n]:
                        dlp_min[n] = dlp2
                        arr1[n][0] = a
                        arr1[n][1] = b
                        arr1[n][2] = c
                        break

    if (count > ncell):
        filename = ('new-sqscell_' + str(ns+1))
        f3 = open('%s.out' % filename, "w")
        f3.write(str(ncell) + '\n\n')
        for i in range(ncell):
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][0][0], arr1[i][0][1], arr1[i][0][2]))
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][1][0], arr1[i][1][1], arr1[i][1][2]))
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][2][0], arr1[i][2][1], arr1[i][2][2])+'\n')
    elif (count > 0):
        filename = ('new-sqscell_' + str(ns+1))
        f3 = open('%s.out' % filename, "w")
        f3.write(str(count) + '\n\n')
        for i in range(count):
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][0][0], arr1[i][0][1], arr1[i][0][2]))
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][1][0], arr1[i][1][1], arr1[i][1][2]))
            f3.write("{:.6f} {:.6f} {:.6f}\n".format(arr1[i][2][0], arr1[i][2][1], arr1[i][2][2])+'\n')
    else:
        print("No cells within the given parameters, change angle or tolerance for %s number of atoms" % str((ns+1)*nat))

if __name__== "__main__":
  main()