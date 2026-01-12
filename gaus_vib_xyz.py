#!/usr/bin/env python
#
# vib-xyz.py v. 0.0.0
#
# Raman off-resonant activity calculator
# using VASP as a back-end.
#
# Contributor: Chloe Groome (UC Irvine)
#
# URL: http://raman-sc.github.io
#
# MIT license, 2013 - 2016
#

import re
import sys, os
from math import sqrt

def parse_gaus(gaus_fh):
    gaus_fh.seek(0) # just in case

    b = []
    symbols = []
    atomtot = 0
    geo = False

    #collect atomic coordinates
    for line in gaus_fh:
        if line == ' \n' and b != []:
            geo = False
            break

        if geo == True:
            atomtot += 1
            symbols.append([x for x in line.split()[0]])
            b.append([float(x) for x in line.split()[1:4]])

        if "Multiplicity" in line:
            geo = True

    
    #collect frequencies 
    eigvals = []

    for line in gaus_fh:
            
        if "Frequencies --" in line:
            print(line.split())
            c = line.split()[2:5]

            eigvals.append(c[0])
            eigvals.append(c[1])
            eigvals.append(c[2])

    #collect eigenvectors
    eigvec1 = [ 0.0 for i in range(atomtot) ]
    eigvec2 = [ 0.0 for i in range(atomtot) ]
    eigvec3 = [ 0.0 for i in range(atomtot) ]

    eigvecs = [ 0.0 for i in range(len(eigvals)) ]
    norms   = [ 0.0 for i in range(len(eigvals)) ]

    eigv = False
    j = 0
    n = 0
    
    #start the line search from the top of the file
    gaus_fh.seek(0)

    for line in gaus_fh:
        c = line.split()[2:]

        if len(c) < 5 and j > atomtot:
            eigv = False
            break

        if eigv == True:
            
            #gaussian outputs 3 eigenvectors per line
            eigvec1[j] = [float(c[0]), float(c[1]), float(c[2])]
            eigvec2[j] = [float(c[3]), float(c[4]), float(c[5])]
            eigvec3[j] = [float(c[6]), float(c[7]), float(c[8])]

            if j == atomtot - 1:
                eigvecs[n]   = eigvec1
                eigvecs[n+1] = eigvec2
                eigvecs[n+2] = eigvec3

                norms[n]   = sqrt( sum([abs(x)**2 for sublist in eigvec1 for x in sublist] ) )
                norms[n+1] = sqrt( sum([abs(x)**2 for sublist in eigvec2 for x in sublist] ) )
                norms[n+2] = sqrt( sum([abs(x)**2 for sublist in eigvec3 for x in sublist] ) )

                n = n + 3
                j = 0
                eigv = False

            else: 
                j = j + 1

        if "  Atom  AN      X      Y      Z" in line:
            eigv = True
    
    #collect dipole deriveratives dQ/dx, dQ/dy, dQ/dz
    irx = []
    iry = []
    irz = []

    #collect raman polarizabilities 
    epsilon = []
    pol = [ 0.0 for i in range(len(eigvals)) ]
    j = 0


    #start the line search from the top of the file
    gaus_fh.seek(0)

    for line in gaus_fh:
            
        if "Dipole derivative wrt mode" in line:
            print(line.split())
            c = line.split()[5:8]

            #convert gaussian string output to python recognizable scientific notation
            irx.append(float(c[0].replace('D', 'e')))
            iry.append(float(c[1].replace('D', 'e')))
            irz.append(float(c[2].replace('D', 'e')))

        if "Polarizability derivatives wrt mode" in line:
            gaus_fh.readline()
            epsilon.append([float(x) for x in line.split()[1:4]])
            epsilon.append([float(x) for x in gaus_fh.readline().split()[1:4]])
            epsilon.append([float(x) for x in gaus_fh.readline().split()[1:4]])

            pol[j] = epsilon
            j += 1


    return b, atomtot, eigvals, eigvecs, norms, irx, iry, irz, pol

# def get_epsilon_from_slurm(slurm_fh):
#     epsilon = []
#     polar1 = [ 0.0 for i in range(len(eigvals)) ]
#     polar2 = [ 0.0 for i in range(len(eigvals)) ]
#     count = 0
#     count1 = 0
#     count2 = 0
#     eps = False

#     slurm_fh.seek(0) # just in case
#     for line in slurm_fh:
        
#         if count == 2*len(eigvals):
#             break

#         if " -----------------------------------------------" in line:
#             eps = False

#         if eps == True:
#             epsilon.append([float(x) for x in line.split()[1:4]])
#             epsilon.append([float(x) for x in slurm_fh.readline().split()[1:4]])
#             epsilon.append([float(x) for x in slurm_fh.readline().split()[1:4]])

#             if count % 2 == 0:
#                 polar2[count2] = epsilon
#                 count2 += 1

#             if count % 2 == 1:
#                 polar1[count1] = epsilon
#                 count1 += 1

#             count += 1
#             epsilon = []

#         if "              X              Y              Z" in line:
#             eps = True
#             slurm_fh.readline() # -----------------------------------------------

#     return polar1, polar2









def parse_env_params(params):
    tmp = params.strip().split('_')
    if len(tmp) != 4:
        print("[parse_env_params]: ERROR there should be exactly four parameters")
        sys.exit(1)
    #
    [first, last, nderiv, step_size] = [int(tmp[0]), int(tmp[1]), int(tmp[2]), float(tmp[3])]
    #
    return first, last, nderiv, step_size

if __name__ == '__main__':
    from math import pi
    from shutil import move
    import os
    import datetime
    import time
    #import argparse
    import optparse
    #
    print("")
    print("    Raman off-resonant activity calculator,")
    print("    using VASP as a back-end.")
    print("")
    print("    Contributors: Alexandr Fonari  (Georgia Tech)")
    print("                  Shannon Stauffer (UT Austin)")
    print("    MIT License, 2013")
    print("    URL: http://raman-sc.github.io")
    print("    Started at: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print("")
    #
    description  = "Before run, set environment variables:\n"
    description += "    VASP_RAMAN_RUN='mpirun vasp'\n"
    description += "    VASP_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]'\n\n"
    description += "bash one-liner is:\n"
    description += "VASP_RAMAN_RUN='mpirun vasp' VASP_RAMAN_PARAMS='1_2_2_0.01' python vasp_raman.py"
    #
    parser = optparse.OptionParser(description=description)
    (options, args) = parser.parse_args()
    #args = vars(parser.parse_args())
    args = vars(options)
    #
    #VASP_RAMAN_RUN = os.environ.get('VASP_RAMAN_RUN')
    VASP_RAMAN_RUN = 'ibrun  vasp_std > out'
    if VASP_RAMAN_RUN == None:
        print("[__main__]: ERROR Set environment variable 'VASP_RAMAN_RUN'")
        print("")
        parser.print_help()
        sys.exit(1)
    print("[__main__]: VASP_RAMAN_RUN='"+VASP_RAMAN_RUN+"'")
    #
    #VASP_RAMAN_PARAMS = os.environ.get('VASP_RAMAN_PARAMS')
    VASP_RAMAN_PARAMS='01_33_2_0.01'

    if VASP_RAMAN_PARAMS == None:
        print("[__main__]: ERROR Set environment variable 'VASP_RAMAN_PARAMS'")
        print("")
        parser.print_help()
        sys.exit(1)
    print("[__main__]: VASP_RAMAN_PARAMS='"+VASP_RAMAN_PARAMS+"'")
    #
    first, last, nderiv, step_size = parse_env_params(VASP_RAMAN_PARAMS)
    assert first >= 1,    '[__main__]: First mode should be equal or larger than 1'
    assert last >= first, '[__main__]: Last mode should be equal or larger than first mode'
    #if args['gen']: assert last == first, "[__main__]: '-gen' mode -> only generation for the one mode makes sense"
    assert nderiv == 2,   '[__main__]: At this time, nderiv = 2 is the only supported'
    disps = [-1, 1]      # hardcoded for
    coeffs = [-0.5, 0.5] # three point stencil (nderiv=2)
    #
    try:
        gaus_fh = open('gaus.out.log', 'r')
    except IOError:
        print("[__main__]: ERROR Couldn't open input file gaus.out.log, exiting...\n")
        sys.exit(1)

    b, atomtot, eigvals, eigvecs, norms, irx, iry, irz, pol = parse_gaus(gaus_fh)
    
    gaus_fh.close()

    output_ir = open('gaus_ir.csv', 'w')
    output_ir.write("# mode,    freq(cm-1),    activity,    activity_x,    activity_y,    activity_z\n")

    output_irx = open('gaus_irx.txt', 'w')
    output_iry = open('gaus_iry.txt', 'w')
    output_irz = open('gaus_irz.txt', 'w')
    
    for i in range(len(eigvals)):
        eigval = float(eigvals[i])
        
        print("")
        print("[__main__]: Mode #%i: frequency %10.7f cm-1;" % ( i+1, eigval))
        
        ir = [0.0 for x in range(3)]


        ir = float(irx[i])**2 + float(iry[i])**2 + float(irz[i])**2
        
        print("")
        print("! %4i,  freq: %10.5f, activity: %10.7f,  activity_x: %10.7f,  activity_y: %10.7f,  activity_z: %10.7f " % (i+1, eigval, ir, irx[i]**2, iry[i]**2, irz[i]**2))
        output_ir.write("%03i,  %10.5f,  %10.7f,  %10.7f,  %10.7f,  %10.7f\n" % (i+1, eigval, ir, irx[i]**2, iry[i]**2, irz[i]**2))
        output_ir.flush()

        output_irx.write("%10.5f  %10.7f\n" % (eigval, irx[i]**2))
        output_irx.flush()

        output_iry.write("%10.5f  %10.7f\n" % (eigval, iry[i]**2))
        output_iry.flush()

        output_irz.write("%10.5f  %10.7f\n" % (eigval, irz[i]**2))
        output_irz.flush()
    output_ir.close()
    output_irx.close()
    output_iry.close()
    output_irz.close()

    bohr2ang = 0.52917724924
    
    output_fh = open('nw_raman.dat', 'w')
    output_fh.write("# mode    freq(cm-1)    alpha    beta2    activity    activity_x    activity_y    activity_z\n")
    for i in range(len(eigvals)):
        eigval = eigvals[i]
        eigvec = eigvecs[i]
        norm = norms[i]
    #     pol1 = polar1[i]
    #     pol2 = polar2[i]
    #     #
        print("")
        print("[__main__]: Mode #%i: frequency %10.7f cm-1; norm: %10.7f" % ( i+1, eigval, norm ))
        
        ra = [[0.0 for x in range(3)] for y in range(3)]
        
    #     for m in range(3):
    #         for n in range(3):
    #             ra[m][n]   = (pol2[m][n] - pol1[m][n]) * (bohr2ang * bohr2ang) * (norm)/(2*step_size)
        
    #             # ra[m][n]   += eps[m][n] * coeffs[j]/step_size * norm * vol/(4.0*pi)
    #         #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3 
        
    #     #
    #     alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
    #     beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0

    #     #x component
    #     alpha_x = (ra[0][0])/3.0
    #     beta2_x = ( (ra[0][0])**2 + (ra[0][0])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2) )/2.0

    #     #y component
    #     alpha_y = (ra[1][1])/3.0
    #     beta2_y = ( (ra[1][1])**2 + (ra[1][1])**2 + 6.0 * (ra[0][1]**2 + ra[1][2]**2) )/2.0

    #     #z component
    #     alpha_z = (ra[2][2])/3.0
    #     beta2_z = ( (ra[2][2])**2 + (ra[2][2])**2 + 6.0 * (ra[0][2]**2 + ra[1][2]**2) )/2.0

    #     print("")
    #     print("! %4i  freq: %10.5f  alpha: %10.7f  beta2: %10.7f  activity: %10.7f  activity_x: %10.7f  activity_y: %10.7f  activity_z: %10.7f " % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2, 45.0*alpha_x**2 + 7.0*beta2_x, 45.0*alpha_y**2 + 7.0*beta2_y, 45.0*alpha_z**2 + 7.0*beta2_z))
    #     output_fh.write("%03i  %10.5f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f\n" % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2, 45.0*alpha_x**2 + 7.0*beta2_x, 45.0*alpha_y**2 + 7.0*beta2_y, 45.0*alpha_z**2 + 7.0*beta2_z))
    #     output_fh.flush()
    #
    # output_fh.close()
