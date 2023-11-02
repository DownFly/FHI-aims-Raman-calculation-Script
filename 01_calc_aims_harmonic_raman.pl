#!/usr/bin/perl
#
#  script to calculate the atomic vibration frequencies for a given atomic structure
#
#  The output will be written as a  .xyz file, which can be
#  interpreted by most common viewing programs (such as jmol, vmd, ... ).
#
#  R. Gehrke/ F. Hanke (2007) 
#
#  Derivative of dipole moment added to get infrared intensities
# 
#  R. Gehrke (2007)
#
#  Derivative of polarizability added to get Raman intensities
#
#  N.Raimbault. (2017)
#
############################################################################################
#
# Begin user-adjusted variables
#
# THESE FOUR VARIABLES NEED CHECKING:
# the AIMS bin directory 

$AIMS_BINDIR = "/thfs1/home/shanghh/gaoyingxiang/05_aims/01_single_node/bin" ;

# aims executable 

$EXE = "aims.191127.scalapack.mpi.x" ;

# calling sequence: any prefixes for serial/parallel execution.
# serial example is default

$CALLING_PREFIX = "yhrun -p thcp1 -N 1 -n 64";
$CALLING_SUFFIX = "";  

#
# example for mpirun with two cores
#
#$CALLING_PREFIX = "mpirun -np 2 ";
#$CALLING_SUFFIX = "";
#
# example for poe with 16 cores:
#
#$CALLING_PREFIX = "poe ";
#$CALLING_SUFFIX = "-procs 16";
#
# End user-adjusted variables
# 
############################################################################################

$HESSIAN_DIAGONALIZATION_EXE    = "vibrations_diagonalizer.191127.scalapack.mpi.x" ;
$EXE_CALL                       = "$CALLING_PREFIX $AIMS_BINDIR/$EXE $CALLING_SUFFIX" ;
$HESSIAN_DIAGONALIZATION_CALL   = "$CALLING_PREFIX $AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE $CALLING_SUFFIX";

# JOBNAME - is the input parameter
$JOBNAME = @ARGV[0];
# delta for finite difference studies, in Angstrom, take it from the second input
$delta = @ARGV[1];
# set reasonable default value for delta
$delta_default = 0.0025;
# define if polarizability is going to be calculated
$polarizability = @ARGV[2];
$polar_default=0; #boolean
$n_atoms_start= @ARGV[3];
$n_atoms_end= @ARGV[4];

print "This is the finite difference vibration calculation for FHI-AIMS\n           Version 190618\n";

if (@ARGV[0] eq '')
{
   $JOBNAME = "basic" ;
}

if (! -e "control.in.$JOBNAME" )
{
   if (-e "control.in")
   {
       system "cp -Rp control.in control.in.$JOBNAME" ;
   }
   else
   {
       print "No suitable control.in file found.\n" ;
       die ;
   }
}

if (! -e "geometry.in.$JOBNAME" )
{
   if (-e "geometry.in")
   {
       system "cp -Rp geometry.in geometry.in.$JOBNAME" ;
   }
   else
   {
       print "No suitable geometry.in file found.\n" ;
       die ;
   }
}

if (! -e "$AIMS_BINDIR/$EXE")
{
    print "Can't find AIMS executable $AIMS_BINDIR/$EXE. \n Aborting. \n";
    die ;
}

if (! -e "$AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE")
{
    print "Can't find Hessian diagonalizer $AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE. Aborting. \n";
    die ;
}


# threshold for taking into account elements in the hessian - negative numbers use all  
# might be required for numerical precision studies 
$HESSIAN_THRESH = -1;

# ... and check if the user specified something
if ($delta == "")
{
    $delta = $delta_default;
    print STDOUT "Using default value for finite difference step size, delta = ",$delta,"!\n";
}
else
{
    print STDOUT "Using user-specified value for finite difference step size, delta = ",$delta,"!\n";
}

if ($polarizability eq '')
{
    $polar = 0; # set all the elements to 0
    print STDOUT "The polarizability will not be calculated. \n";
}
else
{
    $polar_default = 1;
    print STDOUT "The polarizability will be calculated with DFPT. \n";
    print STDOUT "The harmonic Raman spectrum is obtained as explained in Neugebauer et al., J Comput Chem 23: 895–910, 2002 (see Eq.(42))\n";
}

    print STDOUT "\n";
    print STDOUT "!!! IMPORTANT: \n";
    print STDOUT "    Do not simply assume that the default value for delta is good enough.\n";
    print STDOUT "    The balance between too anharmonic (= too large delta) and too much noise\n";
    print STDOUT "    (= too small delta) can be very system-dependent - be sure to check.\n";
    print STDOUT "    delta can be set at the command line for the vibrations script in the following way:";
    print STDOUT "    \n";
    print STDOUT "    <your_path>/bin/aims.vibrations.<version>.pl <your_chosen_name_for_run> <delta> \n";
    print STDOUT "    \n";
    print STDOUT "    for example, for delta=0.002: \n";
    print STDOUT "    \n";
    print STDOUT "    ~/codes/fhi-aims/bin/aims.vibrations.060111.mpi.pl test_0.002 0.002 \n";
    print STDOUT "\n";

# these are the geometry and control input templates
#  
#   fixme: read a job name from parameter list and create the names as 'geometry.in.$(jobname)' etc ... .
#
$INPUTGEOMETRY = join '','geometry.in.',$JOBNAME;
$INPUTCONTROL  = join '','control.in.',$JOBNAME;
$ATOMICMASSES  = join '','masses.',$JOBNAME,'.dat'; 
$XYZFILE       = join '',$JOBNAME,'.xyz';
$HESSIANNAME   = join '','hessian.',$JOBNAME,'.dat'; 
$GRADDIPOLENAME   = join '','grad_dipole.',$JOBNAME,'.dat'; 
$GRADPOLARNAME   = join '','grad_polar.',$JOBNAME,'.dat'; 
$IRNAME   = join '','ir.',$JOBNAME,'.dat'; 
$RAMANNAME   = join '','raman.',$JOBNAME,'.dat'; 

#  This is the output filename prefix&identifier&suffix
$out_prefix = join '',$JOBNAME,'.' ;
$out_suffix = ".out" ;

# internal specifications ... 
$n_atoms    = 0  ;
@species    = () ;
$speciesnum = 0  ;
$mass_number = 0 ;
@species_masses  = () ;
$atommass   = () ;
@coords     = () ;
@energy     = () ;
@constraint = () ;
@lattice_vector = () ;
$periodic = 0 ;

$STARTIN    = join '','startindex.',$JOBNAME,'.dat';
$STARTOUT   = join '','>startindex.',$JOBNAME,'.dat';
$counter    = 0 ; 
$startindex = 0 ;

#  Work begins here 

# check, if there is a startindex (if some jobs were already done) and if so read it
if (-e $STARTIN)
{
    open (INDEX, $STARTIN) ;
    $startindex = <INDEX> ;
    close (INDEX) ;
} 

#  first read geometry-file to get coordinates and do a single calculation
open (INPUT, $INPUTGEOMETRY) ;
$use_constraint = 0 ;
# read one input record ... also works with (<INPUT>)
while ($_ = <INPUT>)
{ 
    # last input fits regular expression ????
    if (/atom/)
    {
        # create an array with all the entries in a line, starting index is 0
	@line = split " ", $_ ; 
        if ( @line[0] eq "atom" )
        {
	    # add coord triple to the end of array coords
	    push @coords, [($line[1], $line[2], $line[3])] ;
	    $species[$n_atoms] = $line[4] ;
	    $n_atoms++ ;
        }
    }
    if (/constraint_region/)
    {
	@line = split " ", $_ ;
	$constraint[$n_atoms-1] = $line[1] ;
	$use_constraint = 1 ;
    }
    if (/lattice_vector/)
    {
        $periodic = 1 ;
        # create an array with all the entries in a line, starting index is 0
	@line = split " ", $_ ; 
        if ( @line[0] eq "lattice_vector" )
        {
	    # create lattice_vector array
	    push @lattice_vector, [($line[1], $line[2], $line[3])] ;
        }
    }
}
close (INPUT) ;


###
# create control-file for displaced structures
# use old density for restart
###

# create control-file
open (INPUT, $INPUTCONTROL) ;
open (OUTPUT, ">control.in") ;
while (<INPUT>)
{
    if ((/relax_geometry/)||(/MD_run/)) 
    {
	# switch off geometry optimization & molecular dynamics
    }
    else 
	{
	    if (/species/) 
	    {
		@line = split " ", $_; 
		$species_name = $line[1];
		print OUTPUT $_;
	    }
	    else 
	    {
		if (/mass/)
		{
		    @line = split " ", $_;
		    $species_mass = $line[1];
		    print OUTPUT $_ ;
		    $mass_number++  ;   # have mass for one more molecule!
		    push @species_masses, [($species_name, $species_mass)];	
		}
		else
		{
		    print OUTPUT $_ ;
		}
	    }
	}
}
# add some other options that are necessary for the vibration calculation.
# after these, there should not be anything possibly wrong ... 
#print OUTPUT "restart_read_only restart.dat\n" ;
print OUTPUT "compute_forces .true. \n";
print OUTPUT "final_forces_cleaned .true. \n" ;

if ($periodic==0)
{
print OUTPUT "output dipole \n" ;
}

# here add if and build "polarizability" file with only zeros!
if ($polar_default==0)
{
  $polar=0;
}
else
{
  # we do not need this in control.in since we have already have. 
  #if ($periodic==0)
  #{
  #  print OUTPUT "DFPT polarizability \n";
  #}
  #else
  #{
  #  print OUTPUT "DFPT dielectric \n";
  #}
}

close (INPUT) ;
close (OUTPUT) ;


# defines the displacements for each atom and direction
@displacements = ( - $delta, $delta ) ;
@coefficients  = (-1, 1) ;
$c_zero = - 1. / (2 * $delta) ;  # force = -dE/dx we want hessian = d^2E/dx^2 = -dF/dx

# displace atoms
print STDOUT "finite difference calculation based on delta = ", $delta, " Angstrom\n" ;

#----------we make parallel here -----------------------
#for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
for ($i_atom = $n_atoms_start; $i_atom <= $n_atoms_end; $i_atom++)
{
    print STDOUT "\n working on atom ",$i_atom,"\n";
    for ($i_coord = 0; $i_coord < 3; $i_coord++)
    {

	# loop through all displacements, the operator $# denotes the length of the 'displacements' array
	for ($i_displacement = 0; $i_displacement <= $#displacements; $i_displacement++)
	{
	    $identifier = join "","i_atom_",$i_atom,".i_coord_",$i_coord,".displ_",$displacements[$i_displacement] ;
	    $geometryfile = join "","geometry.",$identifier,".in" ;
	    
	    # start creating a new geometry.in for this configuration
	    open (OUTPUT, ">$geometryfile") ;
	    
	    print OUTPUT "# temporary structure-file for finite-difference calculation of forces\n" ;
	    print OUTPUT "# displacement " , $displacements[$i_displacement], " of \# atom ", $i_atom, " direction ", $i_coord, "\n" ;
	    
            if ($periodic==1)
            {
	       for ($i_atom_2 = 0; $i_atom_2 < 3; $i_atom_2++)
               {
                   print OUTPUT "lattice_vector ";
	           for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
                   {
                       print OUTPUT $lattice_vector[$i_atom_2][$i_coord_2], " " ;
                   }
                   print OUTPUT "\n" ;
               }
            }
	    for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2++)
	    {
		print OUTPUT "atom " ;
		for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		{
		    # add displacement to coordinate if it is the right one, otherwise use the original coordinate
		    if (($i_coord eq $i_coord_2) && ($i_atom eq $i_atom_2))
		    {
			print OUTPUT $coords[$i_atom_2][$i_coord_2] + $displacements[$i_displacement], " " ;
		    } 
		    else
		    {
			print OUTPUT $coords[$i_atom_2][$i_coord_2], " " ;
		    }
		}
		print OUTPUT $species[$i_atom_2], "\n" ;
		if ($use_constraint == 1) 
		{
		    print OUTPUT "constraint_region ", $constraint[$i_atom_2], "\n" ;
		}
	    }
	    
	    close (OUTPUT) ;
	    $out_file = join "",$out_prefix,$identifier,$out_suffix ;
	    
	    # jump over jobs that are already done - but only jump over the aims calculation, do everything else
	    if (++$counter > $startindex) 
	    {
                system "mkdir $identifier ";
		system "cp control.in  ./$identifier/" ;
		system "mv $geometryfile ./$identifier/geometry.in" ;
		#system "$EXE_CALL > $out_file  < /dev/null " ;
		#system "cp $out_file  ../ " ;
	    }
	    
	    # write next start index: one new calculation is done, write number to file 'startindex.dat'
	    open (INDEX, $STARTOUT) ;
	    print INDEX $counter ;
	    close (INDEX) ;
	    	    

	} # all displacements for a force component done
	
    }     # i_coord    
}         # i_atom


# clean up all the temporary files after completion.
# # All the data files are kept, which means that the hessian diagonalization 
# # and the intensity calculation could be redone from the results. 
#
print "Cleaning temporary files\n\n";

system "cp control.in.$JOBNAME control.in";
system "cp geometry.in.$JOBNAME geometry.in";

