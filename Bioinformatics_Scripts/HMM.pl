$file = shift;
open(IN, $file);
@data = <IN>; 
chomp(@data); 

$seq1 = $data[1]; 
$seq2 = $data[3]; 

#print "$seq1\n";
#print "$seq2\n";

$lens1 = length($seq1);
$lens2 = length($seq2);
#
#print "len 1: $lens1\n";
#print "len 2: $lens2\n";

$delta = shift; 
$eta = shift; 

#### Transition probabilities ####
$BX = $delta; $BY = $delta; $BM = 1 - 2*$delta; 
$MX = $delta; $MY = $delta; $MM = 1 - 2*$delta; 
$XX = $eta; $XY = 0; $XM = 1 - $eta; 
$YY = $eta; $YX = 0; $YM = 1 - $eta; 

#### Emission probabilities ####
$pgap = 1; $pm = 0.6; $pmm = 1 - $pm;

#### Initialize B ####
$B[0][0] = 1; $B[0][1] = 0; $B[1][0] = 0; 
for (my $i=1; $i<=length($seq2);$i++){
  for(my $j=1;$j<=length($seq1);$j++){
    $B[$i][$j] = 0;       
  }  
}

#for(my $i=0; $i<=length($seq2);$i++){
#  for(my $j=0; $j<=length($seq1);$j++){
#  
#    print($B[$i][$j] , " ");
#  }
#  print("\n");
#}

#### Initialize M ####
$M[0][0] = 0; 
for(my $i=1; $i<=length($seq2); $i++){
  $M[$i][0] = 0;
}
for(my $i=1; $i<=length($seq1); $i++){
  $M[0][$i] = 0;
}

#### Initialize X ####
$X[0][0] = 0; 
$X[1][0] = $pgap*$BX*$B[0][0]; 
for(my $i=2; $i<length($seq2); $i++){
  $X[$i][0] = $pgap*$XX*$X[$i-1][0];
}
for(my $i=1; $i<length($seq1); $i++){
  $X[0][$i] = 0; 
}

#### Initialize Y #### 
$Y[0][0] = 0; 
$Y[0][1] = $pgap*$BY*$B[0][0]; 
for(my $i=1; $i<length($seq2); $i++){
  $Y[$i][0] = 0; 
}
for(my $i=2; $i<length($seq1); $i++){
  $Y[0][$i] = $pgap*$YY*$Y[0][$i-1]; 
}

$T[0][0] = 0;
#### Initialize T ####
for(my $i=1;$i<=length($seq2);$i++){
  $T[$i][0] = "U"; 
}
for(my $j=1; $j<=length($seq1);$j++){
   $T[0][$j] = "L"; 
}


#for(my $i=0;$i<=length($seq2);$i++){
#  for(my $j=0;$j<=length($seq1);$j++){
#    print($T[$i][$j] , " ");
#  }
#  print("\n");
#}

$MorMM = -10; 

for(my $i=1; $i<=length($seq2); $i++){
  for(my $j=1; $j<=length($seq1); $j++){
    $ch1 = substr($seq2, $i-1, 1); 
    $ch2 = substr($seq1, $j-1, 1);
    if($ch1 eq $ch2){
      $MorMM = $pm;
    }
    else{
      $MorMM = $pmm; 
    }
    $M1 = $MM*$M[$i-1][$j-1];
    $M2 = $XM*$X[$i-1][$j-1];
    $M3 = $YM*$Y[$i-1][$j-1];
    $M4 = $BM*$B[$i-1][$j-1];  
    $max = max($M1,$M2,$M3,$M4);
    $M[$i][$j] = $MorMM*$max;
    $Y1 = $MY*$M[$i][$j-1];
    $Y2 = $YY*$Y[$i][$j-1]; 
    $Y3 = $BY*$B[$i][$j-1];
    $max = max($Y1,$Y2,$Y3); 
    $Y[$i][$j] = $pgap*$max; 
    $X1 = $MX*$M[$i-1][$j]; 
    $X2 = $XX*$X[$i-1][$j];
    $X3 = $BX*$B[$i-1][$j];
    $max = max($X1,$X2,$X3); 
    $X[$i][$j] = $pgap*$max;
    if($M[$i][$j] >= $X[$i][$j] && $M[$i][$j] >= $Y[$i][$j]){
      $T[$i][$j] = "D"; 
    }
    elsif($X[$i][$j] >= $M[$i][$j] && $X[$i][$j] >= $Y[$i][$j]){
      $T[$i][$j] = "U";
    }
    elsif($Y[$i][$j] >= $M[$i][$j] && $Y[$i][$j] >= $X[$i][$j]){
      $T[$i][$j] = "L";
    }
  }
}



#for(my $i=0; $i<=length($seq2);$i++){
#  for(my $j=0; $j<=length($seq1);$j++){
#    print($T[$i][$j] , " ");
#  }
#  print("\n");
#}



####Traceback 
$aligned_seq1 = "";
$aligned_seq2 = "";
$i = length($seq2);
$j = length($seq1);
while(!($i == 0 && $j == 0)){
	if($T[$i][$j] eq 'D'){
    $ch1 = substr($seq2, $i-1, 1); 
    $ch2 = substr($seq1, $j-1, 1); 
    $i--; $j--; 
   }
  elsif($T[$i][$j] eq 'U'){
    $ch1 = substr($seq2,$i-1,1);
    $ch2 = '-';
    $i--;
  } 
  elsif($T[$i][$j] eq 'L'){
    $ch1 = '-';
    $ch2 = substr($seq1,$j-1,1);
  $j--;
  }
  $aligned_seq1 = $ch1 . $aligned_seq1; 
  $aligned_seq2 = $ch2 . $aligned_seq2; 
}

#print "$data[0]\n"; 
#print "$aligned_seq2\n";
#print "$data[2]\n"; 
#print "$aligned_seq1\n";

$tempfile = "HMMres.txt";

if(-f $tempfile){
      unlink $tempfile;
    }

open(OUTFILE,">>$tempfile") or die;
print OUTFILE "$data[0]\n"; 
print OUTFILE "$aligned_seq2\n";
print OUTFILE "$data[2]\n"; 
print OUTFILE "$aligned_seq1\n"; 
close(OUTFILE); 

#### Max function, returns the max value from array input #### 
sub max{
  my $max = shift(@_);
  foreach $value (@_){
    $max = $value if $max < $value; 
  }
  return $max;
}