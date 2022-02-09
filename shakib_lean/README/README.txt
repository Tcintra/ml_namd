original code given by Farnaz Shakib, farnaz.a.shakib@njit.edu
Shakib & Hanna, J. Chem. Phys. 141, 044122 (2014)


=========================================================================
PCET_model--collective_coord_solvent:
=========================================================================

(model described in shakib_2014.pdf and email--PCET_and_MQCL_codes.pdf)

I) separate all quantum mechanical (q_p, q_e) and classical (Q) components of the Hamiltonian


II) integrate the quantum mechanical part <Psi|H|Psi> "beforehand", expanding Psi in terms of 25 harmonic oscillators for the proton, and 25 harmonic oscillators for the electron

PCET_model--collective_coord_solvent : contains integrations in separate files for 
1) kinetic energy ('Kinetic')
2) potential energy ('Potential')
3) coupling ('Vpes1' and 'Vpes2') 
4) overlap of the wavefunction ('Overlap')


III)

reorder.f90 : reads the output of each submatrix output in each subdirectory (Kinetic, Potential, Vpes1, Vpes2, Overlap) which are based on the indices of proton and electron harmonic oscillators and writes it to a new file with matrix indices from 1 to 625

IV)

matrix.f90 : writes matrix elements for 625 x 625 matrix to matrix.in.  (matrix.in is read in the main code with classical parts incorporated there.)

V) Full matrix is diagonalized at each time step




=========================================================================
mqcl--concerted_PCET_momentum4
=========================================================================

for single classical collective coordinate representing the solvent/environment
(described in email--PCET_and_MQCL_codes.pdf)

src/adiab_nve.f90 : reads matrix.in and builds the model; calculates adiabatic states, Hellman-Feynman forces, nonadiabatic coupling matrix elements (v.d), and d[alpha,beta]




=========================================================================
mqcl--pcet--bath
=========================================================================

for bath of harmonic oscillators representing the classical solvent/environment
(described email--another_code.pdf and in Shakib_2016_rate.pdf)

src/adiab_nve.f90 : reads matrix.in and builds the model; calculates adiabatic states, Hellman-Feynman forces, nonadiabatic coupling matrix elements (v.d), and d[alpha,beta]

force.f90 : calculates classical forces

runnve.f90 : runs velocity-verlet algorithm

