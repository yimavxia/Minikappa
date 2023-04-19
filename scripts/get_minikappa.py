import numpy as np
from conductivity import class_kappa
kappa = class_kappa()
kappa.get_minikappa_phonopy(mesh_in = [12,12,12],                        
                            sc_mat = np.eye(3)*4,
                            pm_mat = np.eye(3),
                            list_temp = [600.0],
                            name_pcell = "POSCAR-unitcell",
                            name_ifc2nd = "FORCE_CONSTANTS",
                            list_taufactor = [2.0])
