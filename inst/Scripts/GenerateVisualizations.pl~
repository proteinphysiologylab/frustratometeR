#----------Create the files and print the headers if needed--------
open(PML, ">$ARGV[2]/$ARGV[1]\_$ARGV[3].pml");
print PML 
"load $ARGV[1]
hide all
unset dynamic_measures
show cartoon, all
color grey, all
run draw_links.py
";

open(JML, ">$ARGV[2]/$ARGV[1]\_$ARGV[3].jml");
print JML "select protein; cartoons;
connect delete;
spacefill off;
color background white;
";

open(TCL, ">$ARGV[2]/$ARGV[1]\_$ARGV[3].tcl");
#----------------------
#----------------------

##-----Read the aux file ------------
my @AuxFile=qx(cat $ARGV[2]/$ARGV[0]);
chomp @AuxFile;

 my $z=0;
 my $y=$z+1;

#-----Write the files---------
for (my $i=0; $i<@AuxFile; $i++)
{
  my @splitted=split " ", $AuxFile[$i];
  print JML "select $splitted[0]:$splitted[2].CA, $splitted[1]:$splitted[3].CA;\nCONNECT single; CONNECT $splitted[5] ";
 
  print TCL "set sel$splitted[0] [atomselect top \"resid $splitted[0] and name CA and chain $splitted[2]\"]\nset sel$splitted[1] [atomselect top \"resid $splitted[1] and name CA and chain $splitted[3]\"]\n# get the coordinates\nlassign [atomselect$z get {x y z}] pos1\nlassign [atomselect$y get {x y z}] pos2\n# draw a green line between the two atoms\ndraw color $splitted[5]\ndraw line \$pos1 \$pos2 style solid width 2\n";

  if($splitted[4] eq "water-mediated")
  {
    if($splitted[5] eq "green")
    {
      print PML "distance min_frst_wm= (///$splitted[2]/$splitted[0]/CA),(///$splitted[3]/$splitted[1]/CA)\n"; 
    }
    else
    {
      print PML "distance max_frst_wm= (///$splitted[2]/$splitted[0]/CA),(///$splitted[3]/$splitted[1]/CA)\n"; 
    }
    print JML "partial radius 0.1\n";
  }
  else
  {
    print PML "draw_links resi $splitted[0] and name CA and Chain $splitted[2], resi $splitted[1] and name CA and Chain $splitted[3], color=$splitted[5], color2=$splitted[5], radius=0.05, object_name=$splitted[0]:$splitted[1]_$splitted[5]\n";
    print JML "single radius 0.1\n";
  }

$z+=2;
$y+=2;
}

#-----Print tails in needed
print PML "zoom all
hide labels
color red, max_frst_wm
color green, min_frst_wm";

print TCL "
mol modselect 0 top all
mol modstyle 0 top newcartoon
mol modcolor 0 top colorid 15
";

close PML;
close JML;
close TCL;

#-------------------

