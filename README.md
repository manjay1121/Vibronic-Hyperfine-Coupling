Program calculates the rate of hyperfine-mediated intersystem crossing for a given system. Program should be executed with Franck-Condon point hyperfine coupling constants (calculated in ORCA), at minimum. For Herzberg-Teller calculations, one can employ phase tracking as an additional calculation. Here, files containing the Q displacements and derivative hyperfine coupling constants are required.

We have now also included a module to incorporate isotopes. Simply edit the block inside the code; turn on the flag, and tell the program the atom number and spin quantum number.

Note: phase tracking is very expensive, and may require submission to the queue. For generation of files, see supplementary scripts. Program of choice is ORCA.

Verbose usage is as follows: ./test14.py -s [SINGLET_HYPERFINE_TABLE] -t [TRIPLET_HYPERFINE_TABLE] -H --qs [<t|Q|s>] --qt [<s|Q}t>] -I --ps [SINGLET_ORBITALS] --pt [TRIPLET_ORBITALS] -p [0-2]

COMMANDS: -s/-t : The isotropic hyperfine coupling constants. As of Orca/6.0.0, ORCA prints these all tablated and rotated into the same axis. NB. this is the minimum requirement for the program to run. -H : Requests a Herzberg-Teller calculation. This is a boolan flag. --qs/--qt : These are the overlap between initial and final state vibrational modes. These can be calculated using the provided FORTRAN script. -I : Requests a Herzberg-Teller calculation with phase tracking. This is a boolean flag. --ps/--pt : The molecular orbital files for both electronic states. Note the REGEX is written to read molden-type formats. Orca can convert .gwb files to molden. -p : The print level between 0-2. 0 is minimal, showing only results. 2 will print all tabulated matrix element components, as well as phase tracking data if appropriate.
