# Build Input file and Cartesian file of Lammps-ReaxFF

Convert XYZMol Format file  to ReaxFF DAT Format file and Build config file Lammps 

**Must be run with:**

      Usage:
          
            perl BuildInputReaxFF.pl [XYZMol-file]
    
      Example (See Folder ./Example/)
          
            perl BuildInputReaxFF.pl CuC74N8H80O8.xyz


*NOTE: To run BuildInputReaxFF.pl the following files are necessary in the working directory:*

    • XYZMol-file          : XYZMol file.

    • BuildInputReaxFF.pl  : The executable file.

    • ReaxFF file          : Reactive MD-force field file of Lammps.

*NOTE: Running Lammps Software see more (https://lammps.sandia.gov/doc/Run_basics.html)*

After a successful run the program will make several output files named as:

    * NameFileXYZ-reaxff.in  : Config file Lammps.
            
    * NameFileXYZ-reaxff.dat : Coords file Lammps.
            
**If you want change options:**

* Convergence criteria

       my $criteria_lammps = "1e-06";

* Specify the maximum number of steps for minimization 
 
       my $steps_lammps      = "500";

* Run Molecular Dynamics (time fs)

       my $steps_run_lammps  = "50000";

* Force Field ReaxFF Lammps

       my $init_relax        = "ffield.reax.Biomolecules";

* Box system and Periodic Boundary conditions

       my $Box_x = 30;

       my $Box_y = 30;

       my $Box_z = 30;
