use strict;
use List::Util qw(shuffle);
use ParsePDB; 
use File::Basename;

my $File=$ARGV[0];
my $jobsDir=$ARGV[1];
my $scriptsDir=$ARGV[2];

my ($pdbfile, $parentdir, $extension) = fileparse($File, qr/\.[^.]*$/);

print "file: $jobsDir/$File\n";

system ("echo Checking_backbone >> $jobsDir/Checkpoints");

#Obtain PDB Chains
print "#Obtain PDB Chains\n";
my @Chains=qx(awk -F \"\" '{if(\$1\$2\$3\$4==\"ATOM\" && \$0~!/DT/ && \$0~!/DG/ && \$0~!/DA/ && \$0~!/DC/) print \$22}' $jobsDir/$File | sort -u);
chomp @Chains;

#Complete check initiates
print "#Complete check initiates\n";
  my $Tag=qx(perl $scriptsDir/VerifComp.pl $jobsDir/$File);
  print "$File is $Tag\n";

  if($Tag eq "Incompleto")
  {
    system ("echo Backbone incomplete >> $jobsDir/Checkpoints");
    print "Completando    python $scriptsDir/MissingAtoms.py $jobsDir/$File; mv $jobsDir/$File\_completed $jobsDir/$File\n";
    system("python $scriptsDir/MissingAtoms.py $jobsDir/$File; mv $jobsDir/$File\_completed $jobsDir/$File");
    system("awk '{if(\$0!~/EXPDTA/) print}' $jobsDir/$File > $jobsDir/aux; mv $jobsDir/aux $jobsDir/$File");
    system ("echo structure completed >> $jobsDir/Checkpoints");
   }
   else
   {
    system ("echo Structure complete >> $jobsDir/Checkpoints");
   }

