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

$CALLING_PREFIX = "yhrun -p thcp1 -N 1 -n 1 ";
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


# start writing MKL output with the original structure 
# also look up the masses of all elements and add them to file $ATOMICMASSES - for input into the diagonalization code
open(MASSFILE,">$ATOMICMASSES");
print "n_atoms = ", $n_atoms, "\n";

for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
{
    # first check through the masses read from control.in whether or not there is one that the script might recognize
    $found = 0;
    for($imasses = 0; $imasses < $mass_number; $imasses++)
    {
	if ($species[$i_atom]=~$species_masses[$imasses][0])
	{
	    $found = 1;
	    $atommass = $species_masses[$imasses][1];
	}
    }    
    if (!$found) 
    {
	# change species into atomic numbers and remember the masses
	# data from http://www.chem.qmul.ac.uk/iupac/AtWt/
	# first, try to match one-letter elements, but some of these might be wrong .... 
	if ($species[$i_atom] =~ /B/  ) { $atommass =  10.811      ;}
	if ($species[$i_atom] =~ /C/  ) { $atommass =  12.0107     ;}
	if ($species[$i_atom] =~ /F/  ) { $atommass =  18.9984032  ;}
	if ($species[$i_atom] =~ /H/  ) { $atommass =  1.00794     ;}
        if ($species[$i_atom] =~ /D/  ) { $atommass =  2.014102    ;} # deuterium
        if ($species[$i_atom] =~ /T/  ) { $atommass =  3.01605     ;} # tritium 
	if ($species[$i_atom] =~ /I/  ) { $atommass =  126.90447   ;}
	if ($species[$i_atom] =~ /K/  ) { $atommass =  39.0983     ;}
	if ($species[$i_atom] =~ /N/  ) { $atommass =  14.0067     ;}
	if ($species[$i_atom] =~ /O/  ) { $atommass =  15.9994     ;}
	if ($species[$i_atom] =~ /P/  ) { $atommass =  30.973762   ;}
	if ($species[$i_atom] =~ /S/  ) { $atommass =  32.065      ;}
	if ($species[$i_atom] =~ /U/  ) { $atommass =  238.02891   ;}
	if ($species[$i_atom] =~ /V/  ) { $atommass =  50.9415     ;}
	if ($species[$i_atom] =~ /W/  ) { $atommass =  183.84      ;}
	if ($species[$i_atom] =~ /Y/  ) { $atommass =  88.90585    ;}
	# ... any errors made so far are fixed by matching two- or three letter elements later
	if ($species[$i_atom] =~ /He/ ) { $atommass =  4.002602    ;}
	if ($species[$i_atom] =~ /Li/ ) { $atommass =  6.941       ;}
	if ($species[$i_atom] =~ /Be/ ) { $atommass =  9.012182    ;}
	if ($species[$i_atom] =~ /Ne/ ) { $atommass =  20.1797     ;}
	if ($species[$i_atom] =~ /Na/ ) { $atommass =  22.98976928 ;}
	if ($species[$i_atom] =~ /Mg/ ) { $atommass =  24.3050     ;}
	if ($species[$i_atom] =~ /Al/ ) { $atommass =  26.9815386  ;}
	if ($species[$i_atom] =~ /Si/ ) { $atommass =  28.0855     ;}
	if ($species[$i_atom] =~ /Cl/ ) { $atommass =  35.453      ;}
	if ($species[$i_atom] =~ /Ar/ ) { $atommass =  39.948      ;}
	if ($species[$i_atom] =~ /Ca/ ) { $atommass =  40.078      ;}
	if ($species[$i_atom] =~ /Sc/ ) { $atommass =  44.955912   ;}
	if ($species[$i_atom] =~ /Ti/ ) { $atommass =  47.867      ;}
	if ($species[$i_atom] =~ /Cr/ ) { $atommass =  51.9961     ;}
	if ($species[$i_atom] =~ /Mn/ ) { $atommass =  54.938045   ;}
	if ($species[$i_atom] =~ /Fe/ ) { $atommass =  55.845      ;}
	if ($species[$i_atom] =~ /Co/ ) { $atommass =  58.933195   ;}
	if ($species[$i_atom] =~ /Ni/ ) { $atommass =  58.6934     ;}
	if ($species[$i_atom] =~ /Cu/ ) { $atommass =  63.546      ;}
	if ($species[$i_atom] =~ /Zn/ ) { $atommass =  65.409      ;}
	if ($species[$i_atom] =~ /Ga/ ) { $atommass =  69.723      ;}
	if ($species[$i_atom] =~ /Ge/ ) { $atommass =  72.64       ;}
	if ($species[$i_atom] =~ /As/ ) { $atommass =  74.92160    ;}
	if ($species[$i_atom] =~ /Se/ ) { $atommass =  78.96       ;}
	if ($species[$i_atom] =~ /Br/ ) { $atommass =  79.904      ;}
	if ($species[$i_atom] =~ /Kr/ ) { $atommass =  83.798      ;}
	if ($species[$i_atom] =~ /Rb/ ) { $atommass =  85.4678     ;}
	if ($species[$i_atom] =~ /Sr/ ) { $atommass =  87.62       ;}
	if ($species[$i_atom] =~ /Zr/ ) { $atommass =  91.224      ;}
	if ($species[$i_atom] =~ /Nb/ ) { $atommass =  92.90638    ;}
	if ($species[$i_atom] =~ /Mo/ ) { $atommass =  95.94       ;}
	if ($species[$i_atom] =~ /Tc/ ) { $atommass =  98          ;}
	if ($species[$i_atom] =~ /Ru/ ) { $atommass =  101.07      ;}
	if ($species[$i_atom] =~ /Rh/ ) { $atommass =  102.90550   ;}
	if ($species[$i_atom] =~ /Pd/ ) { $atommass =  106.42      ;}
	if ($species[$i_atom] =~ /Ag/ ) { $atommass =  107.8682    ;}
	if ($species[$i_atom] =~ /Cd/ ) { $atommass =  112.411     ;}
	if ($species[$i_atom] =~ /In/ ) { $atommass =  114.818     ;}
	if ($species[$i_atom] =~ /Sn/ ) { $atommass =  118.710     ;}
	if ($species[$i_atom] =~ /Sb/ ) { $atommass =  121.760     ;}
	if ($species[$i_atom] =~ /Te/ ) { $atommass =  127.60      ;}
	if ($species[$i_atom] =~ /Xe/ ) { $atommass =  131.293     ;}
	if ($species[$i_atom] =~ /Cs/ ) { $atommass =  132.9054519 ;}
	if ($species[$i_atom] =~ /Ba/ ) { $atommass =  137.327     ;}
	if ($species[$i_atom] =~ /La/ ) { $atommass =  138.90547   ;}
	if ($species[$i_atom] =~ /Ce/ ) { $atommass =  140.116     ;}
	if ($species[$i_atom] =~ /Pr/ ) { $atommass =  140.90765   ;}
	if ($species[$i_atom] =~ /Nd/ ) { $atommass =  144.242     ;}
	if ($species[$i_atom] =~ /Pm/ ) { $atommass =  145         ;}
	if ($species[$i_atom] =~ /Sm/ ) { $atommass =  150.36      ;}
	if ($species[$i_atom] =~ /Eu/ ) { $atommass =  151.964     ;}
	if ($species[$i_atom] =~ /Gd/ ) { $atommass =  157.25      ;}
	if ($species[$i_atom] =~ /Tb/ ) { $atommass =  158.92535   ;}
	if ($species[$i_atom] =~ /Dy/ ) { $atommass =  162.500     ;}
	if ($species[$i_atom] =~ /Ho/ ) { $atommass =  164.93032   ;}
	if ($species[$i_atom] =~ /Er/ ) { $atommass =  167.259     ;}
	if ($species[$i_atom] =~ /Tm/ ) { $atommass =  168.93421   ;}
	if ($species[$i_atom] =~ /Yb/ ) { $atommass =  173.04      ;}
	if ($species[$i_atom] =~ /Lu/ ) { $atommass =  174.967     ;}
	if ($species[$i_atom] =~ /Hf/ ) { $atommass =  178.49      ;}
	if ($species[$i_atom] =~ /Ta/ ) { $atommass =  180.94788   ;}
	if ($species[$i_atom] =~ /Re/ ) { $atommass =  186.207     ;}
	if ($species[$i_atom] =~ /Os/ ) { $atommass =  190.23      ;}
	if ($species[$i_atom] =~ /Ir/ ) { $atommass =  192.217     ;}
	if ($species[$i_atom] =~ /Pt/ ) { $atommass =  195.084     ;}
	if ($species[$i_atom] =~ /Au/ ) { $atommass =  196.966569  ;}
	if ($species[$i_atom] =~ /Hg/ ) { $atommass =  200.59      ;}
	if ($species[$i_atom] =~ /Tl/ ) { $atommass =  204.3833    ;}
	if ($species[$i_atom] =~ /Pb/ ) { $atommass =  207.2       ;}
	if ($species[$i_atom] =~ /Bi/ ) { $atommass =  208.98040   ;}
	if ($species[$i_atom] =~ /Po/ ) { $atommass =  209         ;}
	if ($species[$i_atom] =~ /At/ ) { $atommass =  210         ;}
	if ($species[$i_atom] =~ /Rn/ ) { $atommass =  222         ;}
	if ($species[$i_atom] =~ /Fr/ ) { $atommass =  223         ;}
	if ($species[$i_atom] =~ /Ra/ ) { $atommass =  226         ;}
	if ($species[$i_atom] =~ /Ac/ ) { $atommass =  227         ;}
	if ($species[$i_atom] =~ /Th/ ) { $atommass =  232.03806   ;}
	if ($species[$i_atom] =~ /Pa/ ) { $atommass =  231.03588   ;}
	if ($species[$i_atom] =~ /Np/ ) { $atommass =  237         ;}
	if ($species[$i_atom] =~ /Pu/ ) { $atommass =  244         ;}
	if ($species[$i_atom] =~ /Am/ ) { $atommass =  243         ;}
	if ($species[$i_atom] =~ /Cm/ ) { $atommass =  247         ;}
	if ($species[$i_atom] =~ /Bk/ ) { $atommass =  247         ;}
	if ($species[$i_atom] =~ /Cf/ ) { $atommass =  251         ;}
	if ($species[$i_atom] =~ /Es/ ) { $atommass =  252         ;}
	if ($species[$i_atom] =~ /Fm/ ) { $atommass =  257         ;}
	if ($species[$i_atom] =~ /Md/ ) { $atommass =  258         ;}
	if ($species[$i_atom] =~ /No/ ) { $atommass =  259         ;}
	if ($species[$i_atom] =~ /Lr/ ) { $atommass =  262         ;}
	if ($species[$i_atom] =~ /Rf/ ) { $atommass =  267         ;}
	if ($species[$i_atom] =~ /Db/ ) { $atommass =  268         ;}
	if ($species[$i_atom] =~ /Sg/ ) { $atommass =  271         ;}
	if ($species[$i_atom] =~ /Bh/ ) { $atommass =  272         ;}
	if ($species[$i_atom] =~ /Hs/ ) { $atommass =  270         ;}
	if ($species[$i_atom] =~ /Mt/ ) { $atommass =  276         ;}
	if ($species[$i_atom] =~ /Ds/ ) { $atommass =  281         ;}
	if ($species[$i_atom] =~ /Rg/ ) { $atommass =  280         ;}
	if ($species[$i_atom] =~ /Uub/) { $atommass =  285         ;}
	if ($species[$i_atom] =~ /Uut/) { $atommass =  284         ;}
	if ($species[$i_atom] =~ /Uuq/) { $atommass =  289         ;}
	if ($species[$i_atom] =~ /Uup/) { $atommass =  288         ;}
	if ($species[$i_atom] =~ /Uuh/) { $atommass =  293         ;}
	if ($species[$i_atom] =~ /Uuo/) { $atommass =  294         ;}
    }
    # finally output the mass line - including the geometry and element names for the moment of inertia calculation ... 
    printf MASSFILE "%10.5f %10.4f %10.4f %10.4f %s\n",$atommass,$coords[$i_atom][0],$coords[$i_atom][1],$coords[$i_atom][2],$species[$i_atom];
}
close(MASSFILE);
open (HESSIAN, ">$HESSIANNAME") ;
open (GRAD_DIPOLE, ">$GRADDIPOLENAME") ;
open (GRAD_POLAR, ">$GRADPOLARNAME") ;

# defines the displacements for each atom and direction
@displacements = ( - $delta, $delta ) ;
@coefficients  = (-1, 1) ;
$c_zero = - 1. / (2 * $delta) ;  # force = -dE/dx we want hessian = d^2E/dx^2 = -dF/dx

# displace atoms
print STDOUT "finite difference calculation based on delta = ", $delta, " Angstrom\n" ;

for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
{
    print STDOUT "\n working on atom ",$i_atom,"\n";
    for ($i_coord = 0; $i_coord < 3; $i_coord++)
    {
	# initialise the hessian line for this particular coordinate only
	@hessian = 0. ;

	# initialise derivative of dipole-moment
	@grad_dipole = 0. ;
	@grad_polar = 0. ;

	# loop through all displacements, the operator $# denotes the length of the 'displacements' array
	for ($i_displacement = 0; $i_displacement <= $#displacements; $i_displacement++)
	{
	    $identifier = join "","i_atom_",$i_atom,".i_coord_",$i_coord,".displ_",$displacements[$i_displacement] ;
	    
	    $out_file = join "",$identifier,$out_suffix ;
	    
	    #   now get total energy
	    open (RESULT,"$out_file") ;		
	    $conv = 1. ;
	    $found = 0. ;
	    while (<RESULT>)
	    {
		# extract all forces from output file and make sure that they are converged ... 
		if (/Total\ atomic\ forces/)
		{
		    for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2++)
		    {
			# read line and split it according to the spaces
			$_ = <RESULT> ;
			@line = split " ", $_ ;
			# extract single coordinate forces
			for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
			    {
				$forces[$i_atom_2][$i_coord_2] = @line[$i_coord_2 + 2];
			    }
		    }
		}
		# extract dipole information
		if (/Total\ dipole\ moment/) 
		{
		    # read line and split it according to the spaces
		    @line = split " ", $_ ;
		    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		    {
			$dipole[$i_coord_2] = @line[$i_coord_2 + 6];
		    }
		}
                # extract polarizability information
		if (/Polarizability/) 
		{
		    # read line and split it according to the spaces
		    @line = split " ", $_ ;
		    for ($i_coord_2 = 0; $i_coord_2 < 6; $i_coord_2++)
		    {
			$polar[$i_coord_2] = @line[$i_coord_2 + 2];
		    }
		}

		if (/Have\ a\ nice\ day/)
		{
		    $found = 1 ;
		}
		if (/WARNING\!\ SELF\-CONSISTENCY\ CYCLE\ DID\ NOT\ CONVERGE/)
		{
		    $conv = 0 ;
		}
	    }
	    close (RESULT) ;
	    if (($found eq 0) || ($conv eq 0))
	    {
		print " * WARNING: ",$identifier, " Not converged?\n" ;
		print " Please check this problem before continuing. \n";
		die;
	    }
	    	    
	    #
	    #   calculate and output force derivatives here ... 
	    #
	    #
	    $index = -1;
	    for ($i_atom_2 = 0; $i_atom_2 < $n_atoms;  $i_atom_2++)
	    {
		for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		{
		    $index++;
		    # This is inside the loop over all translations of one single atom, 
		    # add the relevant component of the finite difference to the (already-present) hessian
		    @hessian[$index] += @coefficients[$i_displacement] * $forces[$i_atom_2][$i_coord_2] ;
		} 
	    }
	    
	    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
	    {
		# add the relevant component of the finite difference to the (already-present) gradient of dipole moment
		@grad_dipole[$i_coord_2] += @coefficients[$i_displacement] * $dipole[$i_coord_2] ;
	    }

	    for ($i_coord_2 = 0; $i_coord_2 < 6; $i_coord_2++)
	    {
		# add the relevant component of the finite difference to the (already-present) gradient of polarizability
		@grad_polar[$i_coord_2] += @coefficients[$i_displacement] * $polar[$i_coord_2] ;
	    }

	} # all displacements for a force component done
	
	# calculate correct Hessian line and write it right away
	$index = -1 ;
	for ($i_atom_2 = 0; $i_atom_2 < $n_atoms;  $i_atom_2++)
	{
	    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
	    {
		$index++;
		@hessian[$index] *= $c_zero ; 
		printf HESSIAN "%20.12E ", @hessian[$index];
#		printf STDOUT "%20.12E ", @hessian[$index];
	    }
	} 
	print HESSIAN "\n";
#	print STDOUT "\n";

	for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
	{
	    @grad_dipole[$i_coord_2] *= $c_zero ; 
	    printf GRAD_DIPOLE "%12.5E ", @grad_dipole[$i_coord_2];
#	    printf STDOUT "%12.5E ", @grad_dipole[$i_coord_2];
	}
	print GRAD_DIPOLE "\n";

	for ($i_coord_2 = 0; $i_coord_2 < 6; $i_coord_2++)
	{
	    @grad_polar[$i_coord_2] *= $c_zero ; 
	    printf GRAD_POLAR "%12.5E ", @grad_polar[$i_coord_2];
	}
	print GRAD_POLAR "\n";
    }    
}    
close (HESSIAN) ;
close (GRAD_DIPOLE) ;
close (GRAD_POLAR) ;

# now, call a Fortran diagonalizer and create a final output with the eigenvalues and eigenvectors of the 
# various vibrating modes. 
# necessary parameters: 
#           hesse matrix name
#           grad_dipole name
#           name of atomic mass input file

# need special temporary output from hessian diagonalization == vibrational information in the right format
# $HESSIAN_OUTPUT = join '',$XYZFILE,'.temp';
print   $HESSIAN_DIAGONALIZATION_CALL,' ',$n_atoms,' ',$HESSIANNAME,' ',$GRADDIPOLENAME,' ',$GRADPOLARNAME,' ',$ATOMICMASSES,' ',$HESSIAN_THRESH,' ',$XYZFILE, ' ', $IRNAME, ' ', $RAMANNAME, " | tee ",$JOBNAME,".vib.out ",$HESSIAN_DIAGONALIZATION_SUFFIX,"\n";
system  "$HESSIAN_DIAGONALIZATION_CALL $n_atoms $HESSIANNAME $GRADDIPOLENAME $GRADPOLARNAME $ATOMICMASSES $HESSIAN_THRESH $XYZFILE $IRNAME $RAMANNAME $periodic | tee $JOBNAME.vib.out $HESSIAN_DIAGONALIZATION_SUFFIX";

# clean up all the temporary files after completion.
# All the data files are kept, which means that the hessian diagonalization 
# and the intensity calculation could be redone from the results. 
