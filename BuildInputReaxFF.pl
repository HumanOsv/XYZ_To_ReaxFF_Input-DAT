#!/usr/bin/perl

use strict;
use warnings;

#
# # # 
# Convergence criteria
my $criteria_lammps = "1e-06";
# Specify the maximum number of steps for minimization 
my $steps_lammps      = "500";
# Run Molecular Dynamics (time)
my $steps_run_lammps  = "50000";
# Force Field ReaxFF Lammps
my $init_relax        = "ffield.reax.Biomolecules";
# Box system and Periodic Boundary conditions
my $Box_x = 30;
my $Box_y = 30;
my $Box_z = 30;








# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my $num_atoms_xyz;
##############
# Hashs 
my $other_element = 0.8;
my %Atomic_radii = ( 'H'  => '0.4', 'He' => '0.3', 'Li' => '1.3', 'Be' => '0.9', 
                     'B'  => '0.8', 'C'  => '0.8', 'N'  => '0.8', 'O'  => '0.7', 
					 'F'  => '0.7', 'Ne' => '0.7', 'Na' => '1.5', 'Mg' => '1.3', 
					 'Al' => '1.2', 'Si' => '1.1', 'P'  => '1.1', 'S'  => '1.0', 
					 'Cl' => '1.0', 'Ar' => '1.0', 'K'  => '2.0', 'Ca' => '1.7', 
					 'Sc' => '1.4', 'Ti' => '1.4', 'V'  => '1.3', 'Cr' => '1.3', 
					 'Mn' => '1.4', 'Fe' => '1.3', 'Co' => '1.3', 'Ni' => '1.2', 
					 'Cu' => '1.4', 'Zn' => '1.3', 'Ga' => '1.3', 'Ge' => '1.2', 
					 'As' => '1.2', 'Se' => '1.2', 'Br' => '1.1', 'Kr' => '1.1', 
					 'Rb' => '2.1', 'Sr' => '1.9', 'Y'  => '1.6', 'Zr' => '1.5', 
					 'Nb' => '1.4', 'Au' => '1.4' );

my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',	
	                  '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',  
					  '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',	
                      '107' => 'Bh', '5'   => 'B', 	'35'  => 'Br', '48'  => 'Cd',	
	                  '20'  => 'Ca', '98'  => 'Cf',	'6'   => 'C',  '58'  => 'Ce',	
	                  '55'  => 'Cs', '17'  => 'Cl',	'27'  => 'Co', '29'  => 'Cu',	
	                  '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
	                  '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',	
	                  '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',	
	                  '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',	
	                  '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',	
                      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',	
					  '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
					  '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',	
					  '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',	
                      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',	
					  '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',	
					  '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',	
					  '76'  => 'Os', '8'   => 'O', 	'46'  => 'Pd', '47'  => 'Ag',	
					  '78'  => 'Pt', '82'  => 'Pb',	'94'  => 'Pu', '84'  => 'Po',	
					  '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',	
					  '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',	
					  '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
					  '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
					  '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',	
					  '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',	
					  '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',	
					  '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
					  '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
                      '30'  => 'Zn', '40'  => 'Zr' );

my %Atomic_mass   = ( 'H'   => '1.0079'  ,'He' => '4.003'   ,'Li'  => '6.941'   ,'Be'  => '9.0122',
                      'B'   => '10.811'  ,'C'  => '12.018'  ,'N'   => '14.0067' ,'O'   => '15.9994', 
                      'F'   => '18.998'  ,'Ne' => '20.179'  ,'Na'  => '22.9897' ,'Mg'  => '24.305',
                      'Al'  => '26.981'  ,'Si' => '28.085'  ,'P'   => '30.9738' ,'Cl'  => '35.453',
                      'K'   => '39.098'  ,'Ar' => '39.948'  ,'Ca'  => '40.078'  ,'Sc'  => '44.9559',
                      'Ti'  => '47.867'  ,'V'  => '50.942'  ,'Cr'  => '51.9961' ,'Mn'  => '54.938',
                      'Fe'  => '55.845'  ,'Ni' => '58.693'  ,'Co'  => '58.9332' ,'Cu'  => '63.546',
                      'Zn'  => '65.390'  ,'Ga' => '69.723'  ,'Ge'  => '72.64'   ,'As'  => '74.9216', 
                      'Se'  => '78.960'  ,'Br' => '79.904'  ,'Kr'  => '83.8'    ,'Rb'  => '85.4678', 
                      'Sr'  => '87.620'  ,'Y'  => '88.906'  , 'Zr' => '91.224'  ,'Nb'  => '92.9064',
                      'Mo'  => '95.940'  ,'Tc' => '98.000'  ,'Ru'  => '101.07'  ,'Rh'  => '102.9055',
                      'Pd'  => '106.420' ,'Ag' => '107.868' , 'Cd' => '112.411' ,'In'  => '114.818',
                      'Sn'  => '118.710' ,'Sb' => '121.760' ,'I'   => '126.9045','Te'  => '127.6',
                      'Xe'  => '131.290' ,'Cs' => '132.906' ,'Ba'  => '137.327' ,'La'  => '138.9055',
                      'Ce'  => '140.116' ,'Pr' => '140.908' ,'Nd'  => '144.24'  ,'Pm'  => '145',
                      'Sm'  => '150.360' ,'Eu' => '151.964' ,'Gd'  => '157.25'  ,'Tb'  => '158.9253' ,
                      'Dy'  => '162.500' ,'Ho' => '164.930' , 'Er' => '167.259' ,'Tm'  => '168.9342',
                      'Yb'  => '173.040' ,'Lu' => '174.967' ,'Hf'  => '178.49'  ,'Ta'  => '180.9479',
                      'W'   => '183.840' ,'Re' => '186.207' ,'Os'  => '190.23'  ,'Ir'  => '192.217',
					  'Pt'  => '195.078' ,'Au' => '196.967' ,'Hg'  => '200.59'  ,'Tl'  => '204.3833',
                      'Pb'  => '207.200' ,'Bi' => '208.980' ,'Po'  => '209'     ,'At'  => '210',
					  'Rn'  => '222.000' ,'Fr' => '223.000' ,'Ra'  => '226'     ,'Ac'  => '227',
					  'Pa'  => '231.035' ,'Th' => '232.038' ,'Np'  => '237'     ,'U'   => '238.0289',
					  'Am'  => '243.000' ,'Pu' => '244'     ,'Cm'  => '247'     ,'Bk'  => '247', 
					  'Cf'  => '251.000' ,'Es' => '252'     ,'Fm'  => '257'     ,'Md'  => '258',
					  'No'  => '259.000' ,'Rf' => '261'     ,'Lr'  => '262'     ,'Db'  => '262',
					  'Bh'  => '264.000' ,'Sg' => '266'     ,'Mt'  => '268'     ,'Hs'  => '277' );

					  
###################################
# compute the center of mass
sub measure_center {
	my ($coord_x,$coord_y,$coord_z) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array  = ();
	my $weight = 1;
	# variable sum
	my $sum_weight = 0;
	my $sum_x = 0;
	my $sum_y = 0;
	my $sum_z = 0;
	for ( my $j = 0 ; $j < $num_data ; $j = $j + 1 ){
		$sum_weight+= $weight;
		$sum_x+= $weight * @$coord_x[$j];
		$sum_y+= $weight * @$coord_y[$j];
		$sum_z+= $weight * @$coord_z[$j];		
	}
	my $com_x = $sum_x / $sum_weight;
	my $com_y = $sum_y / $sum_weight;
	my $com_z = $sum_z / $sum_weight;
	# array
	@array = ($com_x,$com_y,$com_z);
	# return array	
	return @array;
}
###################################
# Returns the additive inverse of v(-v)
sub vecinvert {
	my ($center_mass) = @_;
	my @array         = ();
	foreach my $i (@$center_mass) {
		my $invert        = $i * -1;
		$array[++$#array] = $invert; 
	}	
	# return array	
	return @array;
}
###################################
# Returns the vector sum of all the terms.
sub vecadd {
	my ($coord_x,$coord_y,$coord_z,$vecinvert_cm ) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array   = ();
	my $sum_coord_x;
	my $sum_coord_y;
	my $sum_coord_z;
	# array 
	my @array_x = ();
	my @array_y = ();
	my @array_z = ();
	for ( my $i = 0 ; $i < $num_data ; $i = $i + 1 ){	
		$sum_coord_x = @$coord_x[$i]+@$vecinvert_cm[0] ; 
		$sum_coord_y = @$coord_y[$i]+@$vecinvert_cm[1] ;
		$sum_coord_z = @$coord_z[$i]+@$vecinvert_cm[2] ;
		# save array
		$array_x[++$#array_x] = $sum_coord_x;
		$array_y[++$#array_y] = $sum_coord_y;
		$array_z[++$#array_z] = $sum_coord_z;
	}
	@array = ( [@array_x], 
              [@array_y], 
              [@array_z] ); 
	# return array	
	return @array;
}
###################################
# Automatic box size
sub automatic_box_size {
	my ($input_array) = @_;
	my $sum = 0;
	#
	foreach my $i (@{$input_array}) {
		my $radii_val;
		if ( exists $Atomic_radii{$i} ) {
			# exists
			$radii_val = $Atomic_radii{$i};
		} else {
			# not exists
			$radii_val = $other_element ;
		}
		$sum+=$radii_val;
	}
	return $sum;
}
###################################
# read files
sub read_file {
	# filename
	my ($filename) = @_;
	(my $input_file = $filename) =~ s/\s//g;
	my @array          = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	while (my $row = <FILE>) {
		chomp($row);
		push (@array,$row);
	}
	close (FILE);
	# return array	
	return @array;
}
###################################
# delete repeat data
sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}
###################################
# drawing a box around a molecule 
sub box_molecule {
	my ($coordsmin, $coordsmax) = @_;
	#
	my $minx = @$coordsmin[0];
	my $maxx = @$coordsmax[0];
	my $miny = @$coordsmin[1];
	my $maxy = @$coordsmax[1];
	my $minz = @$coordsmin[2];
	my $maxz = @$coordsmax[2];
	# raw the lines
	
	my $filename = 'BOX_kick.vmd';
	open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";	
	print $fh "draw delete all\n";
	print $fh "draw materials off\n";
	print $fh "draw color green\n";
	#
	print $fh "draw line \"$minx $miny $minz\" \"$maxx $miny $minz\" \n";
	print $fh "draw line \"$minx $miny $minz\" \"$minx $maxy $minz\" \n";
	print $fh "draw line \"$minx $miny $minz\" \"$minx $miny $maxz\" \n";
	#
	print $fh "draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\" \n";
	#
	print $fh "draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\" \n";
	#
	print $fh "draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\" \n";
	print $fh "draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\" \n";
	#
	print $fh "draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\" \n";
	print $fh "draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\" \n";
	close $fh;
}
#
sub coords_lammps {
	# Array are send by reference
	my ($frag_1,$num_atoms_xyz,$tmp_box,$init_relax,$file_name) = @_;
	# Reference arrays	
	my $file_xyz = $frag_1;
	my $tmp      = "$file_name-reaxff";
	#
	my $LammpsInput = LammpsInput ($tmp,$file_xyz,"NULL",$init_relax,$tmp_box);
}
###################################
# Config file .in and input file .dat for Molecular Dynamics
sub LammpsInput {
	#
	my $filebase         = $_[0];
	my $coordsMM         = $_[1];
	my $LammpsInput      = "$filebase.in";
	my $LammpsCoords     = "$filebase.dat";
	my $LammpsOutputAxis = "$filebase.lammpstrj";
	my $iteration        = $_[2];
	my $Headerfile       = $_[3];
	my $Box_Length       = $_[4];
	#
	my @words = split (/\n/,$coordsMM);
	my %num_atoms_lammps = ();
	#
	open (COORDSFILE, ">$LammpsCoords");
	print COORDSFILE "# $LammpsCoords file format coords\n";
	print COORDSFILE "\n";
	my $count    = 0;
	my @elements = ();
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		my $label  = $axis[0];  
		push (@elements,$label);
		$count++;
	}
	print COORDSFILE "$count atoms\n";
	my @unique_elements   = uniq @elements;
	my $num_uniq_elements = scalar (@unique_elements);
	print COORDSFILE "$num_uniq_elements atom types\n";
	print COORDSFILE "\n";
	my ($mi_x,$mi_y,$mi_z,$ma_x,$ma_y,$ma_z) = split (" ",$Box_Length);
	print COORDSFILE " $mi_x   $ma_x     xlo xhi\n";
	print COORDSFILE " $mi_y   $ma_y     ylo yhi\n";
	print COORDSFILE " $mi_z   $ma_z     zlo zhi\n";
	print COORDSFILE "\n";
	print COORDSFILE " Masses\n";
	print COORDSFILE "\n";
	my $numb_at = 1;
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		if ( exists $Atomic_mass{$element} ) {
			$mass_val = $Atomic_mass{$element};
		} else {
			$mass_val = $other_element ;
		}
		$num_atoms_lammps{$element} = $numb_at;
		print COORDSFILE " $numb_at $mass_val\n";
		$numb_at++;		
	}
	print COORDSFILE "\n";
	print COORDSFILE " Atoms\n";
	print COORDSFILE "\n";	
	my $count_atoms = 1;
	my $tmp   = "0.0";
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		#
		my $label  = $axis[0];  
		my $axis_x = $axis[1];
		my $axis_y = $axis[2];
		my $axis_z = $axis[3];
		#
		print COORDSFILE "   $count_atoms  $num_atoms_lammps{$label}  $tmp  $axis_x  $axis_y  $axis_z\n";
		$count_atoms++;
	}
	print COORDSFILE "\n";
	close (COORDSFILE);
	#
	open (LAMMPSFILE, ">$LammpsInput");
	print LAMMPSFILE "# REAX FF parameters\n";
	print LAMMPSFILE "#\n";
	print LAMMPSFILE "dimension       3\n";
	print LAMMPSFILE "boundary        p p p\n";
	print LAMMPSFILE "units		    real\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "atom_style      charge\n";
	print LAMMPSFILE "atom_modify     map array sort 0 0.0\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "read_data	    $LammpsCoords\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "pair_style	    reax/c NULL\n";
	#
	print LAMMPSFILE "pair_coeff	    * * $Headerfile ";
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		print LAMMPSFILE "$element ";
	}		
	print LAMMPSFILE "\n";	
	#
	print LAMMPSFILE "neighbor	    2 bin\n";
	print LAMMPSFILE "neigh_modify	every 10 check yes\n";
	print LAMMPSFILE "fix           2 all qeq/reax 1 0.0 10.0 1e-6 reax/c\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "# should equilibrate much longer in practice\n";
	print LAMMPSFILE "#\n";
	print LAMMPSFILE "fix             1 all nvt temp 300 300 10\n";
	print LAMMPSFILE "#fix		        1 all npt temp 273.0 273.0 10.0 iso 1.0 1. 2000.0\n";
	print LAMMPSFILE "timestep        0.2\n";
	print LAMMPSFILE "thermo_style    custom step temp epair pe ke etotal \n";
	print LAMMPSFILE "thermo          1\n";
	print LAMMPSFILE "#\n";    
	print LAMMPSFILE "dump            1 all custom 100 $LammpsOutputAxis id element x y z q \n";
	#
	print LAMMPSFILE "dump_modify      1 element ";
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		print LAMMPSFILE "$element ";
	}		
	print LAMMPSFILE "\n";		
	#
	print LAMMPSFILE "\n";
	print LAMMPSFILE "# $steps_lammps step of minimize\n";	
	print LAMMPSFILE "minimize           $criteria_lammps $criteria_lammps $steps_lammps $steps_lammps\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "# $steps_run_lammps fentoseconds\n";
	print LAMMPSFILE "run                $steps_run_lammps\n";
	#
	return $filebase;
}


################################
# if less than two arguments supplied, display usage
my ($file_name) = @ARGV;
if (not defined $file_name) {
	die "\nConvert Format XYZ to Format ReaxFF DAT; must be run with:\n\nUsage:\n\t perl BuildInputReaxFF.pl [XYZMol-file]\n\n\n";
	exit;  
}
#
my @mol_xyz_file         = read_file ("$file_name");
#
my @mol_cartesian_coords = ();
my @total_atoms          = ();
my @only_coords          = ();
#
for ( my $i=2; $i < scalar(@mol_xyz_file); $i++) {
	my ($atom,@only_coords) = split " ", $mol_xyz_file[$i];
	push (@total_atoms,$atom);
	push (@mol_cartesian_coords,$mol_xyz_file[$i]);
}
#
my $option_box = 1;
#
my @min_coords = ();
my @max_coords = ();
#
if ($option_box == 0) {
	#
	my $sid          = automatic_box_size (\@total_atoms);	
	my $side_box     = sprintf '%.3f',($sid);
	#
	my $side_plus_x  = ($side_box / 2);
	my $side_plus_y  = ($side_box / 2);
	my $side_plus_z  = ($side_box / 2);		
	#
	my $side_minus = (-1 * $side_plus_x);
	@min_coords    = ($side_minus,$side_minus,$side_minus);
	@max_coords    = ($side_plus_x,$side_plus_y,$side_plus_z);
} else {
	my $side_plus_x  = ($Box_x / 2);
	my $side_minus_x = (-1 * $side_plus_x);
	my $side_plus_y  = ($Box_y / 2);
	my $side_minus_y = (-1 * $side_plus_y);
	my $side_plus_z  = ($Box_z / 2);
	my $side_minus_z = (-1 * $side_plus_z);
	#
	@min_coords    = ($side_minus_x,$side_minus_y,$side_minus_z);
	@max_coords    = ($side_plus_x,$side_plus_y,$side_plus_z);
}
#
my $mi_x = sprintf '%.4f',$min_coords[0];
my $mi_y = sprintf '%.4f',$min_coords[1];
my $mi_z = sprintf '%.4f',$min_coords[2];
#
my $ma_x = sprintf '%.4f',$max_coords[0];
my $ma_y = sprintf '%.4f',$max_coords[1];
my $ma_z = sprintf '%.4f',$max_coords[2];
#
print "\n";
print "MESSAGE Box size Min  = $mi_x $mi_y $mi_z\n"; 
print "MESSAGE Box size Max  = $ma_x $ma_y $ma_z\n";
print "MESSAGE ReaxFF Lammps = $init_relax \n";
#
box_molecule (\@min_coords,\@max_coords);
#
(my $without_extension = $file_name) =~ s/\.[^.]+$//;
#
my @coordx = ();
my @coordy = ();
my @coordz = ();
my @elements_new = ();
#
foreach my $i (@mol_cartesian_coords){
	my @Cartesian = split " ", $i;
	push (@coordx,$Cartesian[1]);
	push (@coordy,$Cartesian[2]);
	push (@coordz,$Cartesian[3]);
	#		
	my $element = $Cartesian[0];
	my $radii_val;
	if ( exists $Atomic_number{$element} ) {
		# exists
		$radii_val = $Atomic_number{$element};
	} else {
		# not exists
		$radii_val = $element;
	}
	push (@elements_new,$radii_val);
} 
# call subrutine
my @array_center_mass = measure_center(\@coordx,\@coordy,\@coordz);
my @array_vecinvert   = vecinvert(\@array_center_mass);
# for coords xyz molecules, moveby {x y z} (translate selected atoms)
my @array_catersian   = vecadd (\@coordx,\@coordy,\@coordz,\@array_vecinvert);
#
my $string_coords_total;
for ( my $i = 0 ; $i < scalar (@coordx) ; $i = $i + 1 ){
	my $axis_x = sprintf '%.6f',$array_catersian[0][$i];
	my $axis_y = sprintf '%.6f',$array_catersian[1][$i];
	my $axis_z = sprintf '%.6f',$array_catersian[2][$i];
	#
	$string_coords_total.="$elements_new[$i]  $axis_x  $axis_y  $axis_z\n";
}
$num_atoms_xyz = scalar (@coordx);
my $tmp_box = "$mi_x $mi_y $mi_z $ma_x $ma_y $ma_z";
coords_lammps ($string_coords_total,$num_atoms_xyz,$tmp_box,$init_relax,$without_extension);

print "\n";
print "\n";
print "Files Generated:\n";
print "\n";
print "1 - $without_extension-reaxff.in  (Config file Lammps)\n";
print "2 - $without_extension-reaxff.dat (Coords file Lammps)\n";
print "\n";
print "Run Lammps: \n";
print "         laptop-user\$ lmp_serial -in $without_extension-reaxff.in \n";
print "\n";
print "\n";
print "* * Normal Termination * *\n";
print "\n";
print "\n";
