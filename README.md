# Build Input file and Cartesian file of Lammps

Convert XYZMol Format file  to ReaxFF DAT Format file and Build config file Lammps 

**Must be run with:**

      Usage:
          
            perl BuildInputReaxFF.pl [XYZMol-file]
    
      Example (See Folder ./Example/)
          
            perl BuildInputReaxFF.pl CuC74N8H80O8.xyz
    

After a successful run the program will make several output files named as:

            NameFileXYZ-reaxff.in  (Config file Lammps)
            
            NameFileXYZ-reaxff.dat (Coords file Lammps)
            
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

       my $Box_x = 15;

       my $Box_y = 15;

       my $Box_z = 15;
