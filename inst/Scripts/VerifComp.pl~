use strict;

sub Verif
{

my $completed="true";
my @test=();
my $resi_test=0;
my $missedC=0;
my $missedCA=0;
my $missedCB=0;
my $missedN=0;
my $missedO=0;
my $resAnt="";
my %Gly_res=();

my ($PDB)=@_;


my $iRes=();
my @idRES=();

my $nRES=();
my @NRES=qx(tail -n 2 $PDB);
if(substr($NRES[1], 0,4) eq "ATOM")
{
 $nRES=substr($NRES[1], 23, 3);
}
else
{
 $nRES=substr($NRES[0], 23, 3);
}

$nRES=~s/\s//;


#--------------------------------------
my $nCA=qx(grep \"CA  \\w[a-zA-Z]\" $PDB | wc);
my @CAs=split(" ", $nCA);
$nCA=$CAs[0];


if($nCA ne $nRES)
{
 $completed="false";
}
#--------------------------------------
my $nCB=qx(grep \"CB  \\w[a-zA-Z]\" $PDB |wc);
my @CBs=split(" ", $nCB);
$nCB=$CBs[0];

my $nGLY=qx(grep \"CA  GLY\" $PDB | wc);
my @GLYs=split(" ", $nGLY);
$nGLY=$GLYs[0];

my @Glys=qx(grep \"CA  GLY\" $PDB);

for (my $i=0; $i<@Glys; $i++)
{
	my $Gly_index=substr($Glys[$i], 22, 4);

		$Gly_index =~ s/^\s+//;
		$Gly_index =~ s/\s+$//;
		$Gly_index=$Gly_index;
		if(not exists($Gly_res{$Gly_index})){$Gly_res{$Gly_index}="gly"; }
}

if($nCB ne ($nRES-$nGLY))
{
 $completed="false";
}

#--------------------------------------
my $nC=qx(grep \"C   \\w[a-zA-Z]\" $PDB |wc);
my @Cs=split(" ", $nC);
$nC=$Cs[0];


if($nC ne $nRES)
{
 $completed="false";
}

#--------------------------------------
my $nO=qx(grep \"O   \\w[a-zA-Z]\" $PDB |wc);
my @Os=split(" ", $nO);
$nO=$Os[0];


if($nO ne $nRES)
{
 $completed="false";
}

#--------------------------------------
my $nN=qx(grep \"N   \\w[a-zA-Z]\" $PDB |wc);
my @Ns=split(" ", $nN);
$nN=$Ns[0];

if($nN ne $nRES)
{
 $completed="false";
}

if ($completed eq "true")
{
#print "El pdb esta completo\n";
}
else
{
#print "El pdb Esta incompleto\n";
}
return $completed;
}
#--------------------------------

my $status= &Verif("$ARGV[0]");

if ($status eq "true")
{
print "OK";
}
else
{
print "Incompleto";
}



