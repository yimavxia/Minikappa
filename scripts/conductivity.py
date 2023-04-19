import numpy as np
import math
import cmath
import os
import sys

global debug

class class_kappa:
    """
    input: obj_poscar
    """
    def __init__(self):
        pass

    def get_minikappa_phonopy(self,
                          mesh_in = [8,8,8],
                          sc_mat = np.eye(3)*2,
                          pm_mat = np.eye(3),
                          list_temp = [300.0],
                          name_pcell = "POSCAR-unitcell",
                          name_ifc2nd = "FORCE_CONSTANTS",
                          list_taufactor = [2.0]):
        try:
            import phonopy
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
            from phonopy.interface.calculator import read_crystal_structure
        except:
            print ("Phonopy API version of 2.7.1 is required!")
            
        # load phonopy object
        phonon = phonopy.load(supercell_matrix = sc_mat,
                              primitive_matrix = pm_mat,
                              unitcell_filename = name_pcell,
                              is_symmetry = False,
                              force_constants_filename = name_ifc2nd)
        primcell = (phonon.get_primitive()).cell.T
        volpc = np.abs(np.dot(np.cross(primcell[1],primcell[2]),primcell[0]))
        print (f"Number of atoms in primitive cell: {len(phonon.primitive.masses)}")

        # phonon
        phonon.run_mesh(mesh_in,
                        with_eigenvectors = False,
                        is_gamma_center = True,
                        with_group_velocities = True,
                        is_time_reversal = False,
                        is_mesh_symmetry = False)
        mesh_dict = phonon.get_mesh_dict()
        qpoints = mesh_dict['qpoints']
        weights = mesh_dict['weights']
        freqs = mesh_dict['frequencies']
        eigs = mesh_dict['eigenvectors']
        gvfull = mesh_dict['group_velocities']
        
        # unit
        hbar = 1.054571726470000E-022
        kB = 1.380648813000000E-023
        pi = np.pi
        volpc = volpc/1000.0 
        gvfull = gvfull/10.0 
        freqs = freqs*2*pi   
        nband = len(freqs[0])
        freqcf = 0.1 

        # kappa
        for temp in list_temp:
            for factor in list_taufactor:
                kappaband=np.zeros((nband,nband,3,3), dtype=np.complex128, order='C')
                nqpt=len(qpoints)
                nband=len(freqs[0])
                for i in range(nband):
                    for j in range(nband):
                        for iq in range(nqpt):
                            for k in range(3):
                                for kp in range(3):
                                    omega1=freqs[iq,i]
                                    omega2=freqs[iq,j]
                                    if omega1>freqcf and omega2>freqcf:
                                        if (freqs[iq,i]/2/pi) > 0:
                                            Gamma1=freqs[iq,i]/2/pi*factor
                                        else:
                                            Gamma1=1E10
                                        if (freqs[iq,j]/2/pi) > 0:
                                            Gamma2=freqs[iq,j]/2/pi*factor
                                        else:
                                            Gamma2=1E10
                                        fBE1=1.0/(np.exp(hbar*omega1/kB/temp)-1.0)
                                        fBE2=1.0/(np.exp(hbar*omega2/kB/temp)-1.0)
                                        tmpv=(gvfull[iq,i,j,k]*gvfull[iq,j,i,kp]).real
                                        kappaband[i,j,k,kp]=kappaband[i,j,k,kp]+(omega1+omega2)/2* \
                                            (fBE1*(fBE1+1)*omega1+fBE2*(fBE2+1)*omega2)*tmpv \
                                            /(4*(omega1-omega2)**2+(Gamma1+Gamma2)**2)*\
                                            (Gamma1+Gamma2)
                # conversion
                print (f"Factor: {factor}")
                kappaband=kappaband*1E21*hbar**2/(kB*temp*temp*volpc*nqpt)
                kappaD  = np.zeros((3,3), dtype=np.complex128, order='C')
                kappaOD = np.zeros((3,3), dtype=np.complex128, order='C')
                kappaF  = np.zeros((3,3), dtype=np.complex128, order='C')
                for i in range(nband):
                    for j in range(nband):
                        kappaF=kappaF+kappaband[i,j]
                        if i==j:
                            kappaD=kappaD+kappaband[i,j]
                        else:
                            kappaOD=kappaOD+kappaband[i,j]
                f = open("minikappa"+"-"+str(temp)+"-"+str(factor)+".dat", "w")
                f.write("   ".join(map(str, np.round(kappaD.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.write("   ".join(map(str, np.round(kappaOD.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.write("   ".join(map(str, np.round(kappaF.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.close()
                if True:
                    print ("Diagonal part of minimum thermal conductivity: ")
                    print (kappaD.real[0,0])
                    print ("Off-diagonal part of minimum thermal conductivity: ")
                    print (kappaOD.real[0,0])
                    print ("Total minimum thermal conductivity: ")
                    print (kappaF.real[0,0])
        #return [kappaD.real[0,0], kappaOD.real[0,0]]
