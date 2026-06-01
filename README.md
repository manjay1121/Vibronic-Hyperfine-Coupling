UPDATE: This work is now published [https://doi.org/10.1063/5.0326602]

This program calculates the rate of hyperfine-mediated intersystem crossing for a given system. At a minimum, this program requires hyperfine coupling constants evaluated at the Franck-Condon point. The code is written to interpret formatting exported by ORCA, however this can easily be reconfigured. For Herzberg-Teller calculations, files containing the nuclear displacements (Q) and the derivatives of the hyperfine coupling constants with respect to those Q are required. Additionally, one should employ normal mode phase tracking as discussed in the publication.

Note: phase tracking is very expensive and may require submission to a queue rather than being executed on a local machine. For generation of files, see supplementary scripts. This program was originally built for compatibility with ORCA and the MOLDEN output for orbitals. Orca can convert .gwb files to .molden files. Thankfully this is a prolific format.

We have now also included a module to incorporate a single variation in nuclear spin (to allow for isotopic labelling). This is called as a command line arguement.  

Verbose usage is as follows: ./Project:Superfine.py 
--dE Energy gap
--dEu Splitting in energy gap (upper)
--dEl Splitting in energy gap (lower)
--lst Reorganisation energy from S-T
--lts Reorganisation energy from T-S
--lu Splitting in reorganisation energy (upper)
--ll Splitting in reorganisation energy (lower)
--zp Zeeman pulling effect (if applicable)
-T Temperature
--dH Derivative splitting
-s [SINGLET_HYPERFINE_TABLE] 
-t [TRIPLET_HYPERFINE_TABLE] 
--ht 
--qs [<t|Q|s>] 
--qt [<s|Q}t>] 
--phase 
--ps [SINGLET_ORBITALS] 
--pt [TRIPLET_ORBITALS] 
--iso Request isotope
--index Atom number, as per hyperfine calculation NOT molecular definition
--spin Isotope spin
-p [0-2]

COMMANDS: 
--dE : The adiabatic energy gap between the singlet and triplet biradical spin-states. This should be in units of eV. Sign is important, noted in this implementation as a negative energy corresponding to the triplet being higher in energy than the singlet. 
--dEu/--dEl : Splitting in the adiabatic energy gap due to ZFS and/or MFEs. Should be in units of eV. Is specifically different with respect to the upper/lower shift with respect to the position of the |0> level. 
--lst/--lts : The reorganisation energy between the singlet to triplet state, or the triplet to singlet state, respectively. If no [--lts] input is given the code will assume lts==lst. This value must always be nonzero!. 
--lu/-ll : Splitting in the reorganisation energy due to ZFS and/or MFEs. Should be units of eV. 
--zp : Zeeman pulling effect. If given, this parameter will pull/push the energy of the singlet manifold accordingly. Sign convention dicates the direction of pulling. This parameter is typically very weak to negligible unless in large magnetic fields. 
-T : Temperature of the system, in K. 
--dH The normal mode projection magnitude. How much the derivative vector is projected for each normal mode, with respect to the Franck-Condon point. Is dimensionless; a reasonable value is 0.05. Should be the same as the number used to calculate the hyperfine derivatives.   
-p [0-2]
-s/-t : The isotropic hyperfine coupling constants. As of Orca/6.0.0, ORCA prints these tabulated and rotated into the same axis. N.B. this is the minimum requirement for the program to run. 
--ht (bool): Requests a Herzberg-Teller calculation. 
--qs/--qt : These are the overlap between initial and final state vibrational modes. These can be calculated using the provided FORTRAN script which reads from provided out from the .hess file from ORCA. 
--phase (bool): Requests a Herzberg-Teller calculation with phase tracking. 
--ps/--pt : The molecular orbital files for both electronic states, at the Franck-Condon point. 
--iso (bool) : Requests an isotope substitution. 
--index : The index according to which atom should be substituted for an isotope. This index is given in the basis of the hyperfine calculation (the order in which each atom appears in the hyperfine results), NOT molecular definition
--spin : The new spin of atom [index]. Can be given as a fraction e.g. 5/2. 
-p : The print level (0-2). 0 is minimal, showing only results. 2 will print all tabulated matrix element components as well as phase tracking output.

Provided alongside the code is the [HT] directory. This contains the tabulated data files that the codebase reads and uses, as well as the subcodes required to generate each datapoint if so required. To do so, one requires the hyperfine coupling constants derived from nuclear displacements using each normal mode as a displacement vector; to gain these, proceed into [SINGLET]/[TRIPLET] and run the <normal_mode_propagator_orca> bash code, which will proceed to submit the prerequisite calculations ~(3N-6)*2 per manifold. This will read the frequency calculation [FMN.out], perturb the geometry along each projection using [reader.py], then submit the calculations. Note, you will need to reference your submission script in [normal_mode_propagator_orca]. Following the conclusion of each, run the [hyperfine_extractor] code to pull the Euler-rotated hyperfine coupling constants for [Project:Superfine] to read. Following this, one can run the [hyperfine_cleaner] to purge all non-neccesary files from the directory. Amend the file handes as required. Note for our system, this entire folder came out to >60 GB. Finally, to generate <f|Q|i>, you can use the code provided in [Q_factors]. Written by a dear colleage Dr. Zifei Chen, this bootstrapped Fortran code reads output taken directly from the ORCA .hess files (between initial and final states) to calculate the eigenvector overlap for each normal mode. The variable <natoms> in the code [HR.f90] will need to be amended, then recompiled into [a.out]. Check [README.00]. The resulting file gives a 3N-6 x 5 table which reports the overlaps, the energies of each normal, and the corresponding Huang-Rhys factors. 
