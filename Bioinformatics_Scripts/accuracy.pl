$true_file=shift;
open(IN, $true_file);
@data=<IN>;
chomp(@data);

$t_seq1 = $data[1];
$t_seq2 = $data[3];

#print "t_seq1:\n";
#print "$t_seq1\n";
#print "t_seq2:\n";
#print "$t_seq2\n";


open(IN, shift);
@data =<IN>;

$c_seq1 = $data[1];
$c_seq2 = $data[3];

#print "c_seq1: ";
#print "$c_seq1\n";
#print "c_seq2: ";
#print "$c_seq2\n";


%t_hash = ();
$pos1 = 1;
$pos2 = 1;

for($i=0; $i < length($t_seq1); $i++){
  $ch1 = substr($t_seq1, $i, 1);
  $ch2 = substr($t_seq2, $i, 1);
  if($ch1 ne "-" && $ch2 ne "-" && $ch1 =~ /[A-Z]/ && $ch2 =~ /[A-Z]/){
    $t_hash{$pos1} = $pos2;
  }
  if(substr($t_seq1, $i, 1) ne "-"){
    $pos1++;
  }
  if(substr($t_seq2, $i, 1) ne "-"){
    $pos2++;
  }
}

@t_k = keys(%t_hash); # return keys
$den = scalar(@t_k); # length of list
$pos1 = 1;
$pos2 = 1;

for($i=0; $i < length($c_seq1); $i++){
  if(substr($c_seq1, $i, 1) ne "-" && substr($c_seq2, $i, 1) ne "-"){
    if($t_hash{$pos1} == $pos2){$num++;}
  }
  if(substr($c_seq1, $i, 1) ne "-"){
    $pos1++;
  }
  if(substr($c_seq2, $i, 1) ne "-"){
    $pos2++;
  }
}

$accuracy = ($num/$den);
#print "accuracy: ";
#print "$accuracy\n";

#print "num: ";
#print "$num\n";
#print "den: ";
#print "$den\n";

$tempfile = "acc.txt";

if(-f $tempfile){
      unlink $tempfile;
    }

open(OUTFILE,">>$tempfile") or die;
print OUTFILE "$accuracy\n";  
close(OUTFILE); 
