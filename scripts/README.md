# Python Scripts for Computing Minimum Lattice Thermal Conductivity

Here, we include:

* An implementation of minimum lattice thermal conductivity: conductivity.py
* A modified Phonopy source code to compute off-diagonal group velocities: group_velocity.py
  * To use this script, replace the original Phonopy group_velocity.py with the modified one
* A script to compute minimum lattice thermal conductivity using data in the example folder : get_minikappa.py
  * To run the example, go to `example` and run `python ../get_minikappa.py`
  * You might want to update the input parameters for your specific compounds
  * Only structure file (e.g., POSCAR-unitcell) and harmonic force constants (e.g., FORCE_CONSTANTS) are required to compute KL_min
