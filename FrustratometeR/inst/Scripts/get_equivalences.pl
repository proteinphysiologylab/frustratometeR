#!/usr/local/bin/perl
use strict;


open (SALIDA, ">$ARGV[0]_equivalences.txt");


my %Chains=();
my @Residues= qx(awk '{print}' $ARGV[0]);

my $aux="";

my $eq=1;


#get chains
my $aux_chain='1';
foreach my $line (@Residues)
{
			if(substr($line, 0,4) eq "ATOM")
			{
				my $chain=substr($line, 21,1);
				if ($chain ne $aux_chain)
				{
						$Chains{$chain}=1;
				}
			$aux_chain=$chain;
	}
}

foreach my $single_chain (sort keys %Chains)
{
 my @Residues= qx(awk 'BEGIN{FS=""} { if ((\$1\$2\$3\$4=="ATOM" || (\$1\$2\$3\$4\$5\$6=="HETATM" && \$18\$19\$20=="MSE")) && \$22=="$single_chain" && \$19\$20!~/DT/ && \$19\$20!~/DG/ && \$19\$20!~/DA/ && \$19\$20!~/DC/) print }' $ARGV[0]);

for(my $i=0; $i<@Residues; $i++)
{
			if(substr($Residues[$i], 0,4) eq "ATOM" || (substr($Residues[$i], 17,3) eq "MSE" && substr($Residues[$i], 0,6) eq "HETATM" ) )
			{
				my $res=substr($Residues[$i], 22, 5);
				my $chain=substr($Residues[$i], 21, 1);

				$res =~ s/\s+$//;

				if($res ne $aux)	
				{
					if($eq<10){$eq="   ".$eq;} elsif($eq<100){$eq="  ".$eq;}elsif($eq<1000){$eq=" ".$eq;}
					print SALIDA "$chain $eq $res\n";
					$eq++;
				}
				$aux=$res;
		}
}

}

close SALIDA;
