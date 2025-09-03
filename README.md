This program calculates the rate of hyperfine-mediated intersystem crossing for a given system. At a minimum, this program requires hyperfine coupling constants evaluated at the Franck-Condon point. The code is written to interpret formatting exported by ORCA, however this can easily be reconfigured. For Herzberg-Teller calculations, files containing the nuclear displacements (Q) and the derivatives of the hyperfine coupling constants with respect to those Q are required. Additionally, one should employ normal mode phase tracking as discussed in the publication.

Note: phase tracking is very expensive and may require submission to a queue rather than being executed on a local machine. For generation of files, see supplementary scripts. This program was originally built for compatibility with ORCA and the MOLDEN output for orbitals. Orca can convert .gwb files to .molden files. Thankfully this is a prolific format.

We have now also included a module to incorporate a single variation in nuclear spin (to allow for isotopic labelling). Go to the ISOTOPIC DEFINITIONS block inside the code, turn on the flag, and tell the program the atom number and nuclear spin quantum number of the atom which needs to be labelled. NB. code can be amended to incorporate multiple variations if required. 

Verbose usage is as follows: ./Project:Superfine.py 
-s [SINGLET_HYPERFINE_TABLE] 
-t [TRIPLET_HYPERFINE_TABLE] 
-H 
--qs [<t|Q|s>] 
--qt [<s|Q}t>] 
-I 
--ps [SINGLET_ORBITALS] 
--pt [TRIPLET_ORBITALS] 
-p [0-2]

COMMANDS: 
-s/-t : The isotropic hyperfine coupling constants. As of Orca/6.0.0, ORCA prints these tabulated and rotated into the same axis. N.B. this is the minimum requirement for the program to run. 
-H (bool): Requests a Herzberg-Teller calculation. 
--qs/--qt : These are the overlap between initial and final state vibrational modes. These can be calculated using the provided FORTRAN script which reads from provided out from the .hess file from ORCA. 
-I (bool): Requests a Herzberg-Teller calculation with phase tracking. 
--ps/--pt : The molecular orbital files for both electronic states. 
-p : The print level (0-2). 0 is minimal, showing only results. 2 will print all tabulated matrix element components as well as phase tracking output.

Provided alongside the code is the [HT] directory. This contains the tabulated data files that the codebase reads and uses, as well as the subcodes required to generate each datapoint if so required. To do so, one requires the hyperfine coupling constants derived from nuclear displacements using each normal mode as a displacement vector; to gain these, proceed into [SINGLET]/[TRIPLET] and run the <normal_mode_propagator_orca> bash code, which will proceed to submit the prerequisite calculations ~(3N-6)*2 per manifold. This will read the frequency calculation [FMN.out], perturb the geometry along each projection using [reader.py], then submit the calculations. Note, you will need to reference your submission script in [normal_mode_propagator_orca]. Following the conclusion of each, run the [hyperfine_extractor] code to pull the Euler-rotated hyperfine coupling constants for [Project:Superfine] to read. Following this, one can run the [hyperfine_cleaner] to purge all non-neccesary files from the directory. Amend the file handes as required. Note for our system, this entire folder came out to >60 GB. Finally, to generate <f|Q|i>, you can use the code provided in [Q_factors]. Written by a dear colleage Dr. Zifei Chen, this bootstrapped Fortran code reads output taken directly from the ORCA .hess files (between initial and final states) to calculate the eigenvector overlap for each normal mode. The variable <natoms> in the code [HR.f90] will need to be amended, then recompiled into [a.out]. Check [README.00]. The resulting file gives a 3N-6 x 5 table which reports the overlaps, the energies of each normal, and the corresponding Huang-Rhys factors. 
