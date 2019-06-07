# Build Input file and Cartesian file of Lammps


If you want change options:
* Convergence criteria

       my $criteria_lammps = "1e-06";

* Specify the maximum number of steps for minimization 
 
       my $steps_lammps      = "500";

* Run Molecular Dynamics (time)

       my $steps_run_lammps  = "50000";

* Force Field ReaxFF Lammps

       my $init_relax        = "ffield.reax.Biomolecules";

* Box system and Periodic Boundary conditions

       my $Box_x = 15;

       my $Box_y = 15;

       my $Box_z = 15;
