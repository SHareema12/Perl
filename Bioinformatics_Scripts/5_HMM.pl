$file = shift;
open(IN, $file);
@data = <IN>; 
chomp(@data); 

$seq1 = $data[1]; 
$seq2 = $data[3]; 

$delta = .2; 
$eta = .4; 

@B = (); @M = (); @X = (); @Y = (); @f = (); @b = (); @p = ();

#### Emission probabilities ####
$pgap = 1; $pm = 0.6; $pmm = 1 - $pm; $tau = 0; 

#### Transition probabilities ####
$BX = $delta; $BY = $delta; $BM = 1 - 2*$delta; 
$MX = $delta; $MY = $delta; $MM = 1 - 2*$delta - $tau; 
$XX = $eta; $XY = 0; $XM = 1 - $eta - $tau; 
$YY = $eta; $YX = 0; $YM = 1 - $eta - $tau;

#### Initialization ####
$B[0][0] = 1; $B[0][1] = 0; $B[1][0] = 0;
$M[0][0] = 0; $M[0][1] = 0; $M[1][0] = 0;
$X[1][0] = $pgap*$BX*$B[0][0];
$Y[0][1] = $pgap*$BY*$B[0][0];
$X[0][0] = 0; $X[0][1] = 0;
$Y[0][0] = 0; $Y[1][0] = 0; 
$T[0][1] = 'L'; $T[1][0] = 'U'; 

#for(my $i=2; $i<=length($seq2); $i++){
#  $M[$i][0] = 0; $X[$i][0] = $pgap*$XX*X[$i-1][0]; 
#  $Y[$i][0] = 0; $T[$i][0] = 'U';  
#}
#for(my $j=2; $j<=length($seq1); $j++){
#  $M[0][$j] = 0; $X[0][$j] = 0; 
#  $Y[0][$j] = $pgap*$YY*Y[0][$j-1]; $T[0][$j] = 'L';  
#}

#### Initialize M ####
for(my $i=1; $i<=length($seq2); $i++){
  $M[$i][0] = 0;
}
for(my $i=1; $i<=length($seq1); $i++){
  $M[0][$i] = 0;
}

#### Initialize X ####
$X[1][0] = $pgap*$BX*$B[0][0]; 
for(my $i=2; $i<length($seq2); $i++){
  $X[$i][0] = $pgap*$XX*$X[$i-1][0];
}
for(my $i=1; $i<length($seq1); $i++){
  $X[0][$i] = 0; 
}

#### Initialize Y #### 
$Y[0][1] = $pgap*$BY*$B[0][0]; 
for(my $i=1; $i<length($seq2); $i++){
  $Y[$i][0] = 0; 
}
for(my $i=2; $i<length($seq1); $i++){
  $Y[0][$i] = $pgap*$YY*$Y[0][$i-1]; 
}

#### Initialize T ####
for(my $i=1;$i<=length($seq2);$i++){
  $T[$i][0] = 'U'; 
}
for(my $j=1; $j<=length($seq1);$j++){
   $T[0][$j] = 'L'; 
}


#### Forward Recurrence ####
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
    $M[$i][$j] = $MorMM*$M1 + $MorMM*$M2 + $MorMM*$M3 + $MorMM*$M4;
    $Y1 = $MY*$M[$i][$j-1];
    $Y2 = $YY*$Y[$i][$j-1]; 
    $Y3 = $BY*$B[$i][$j-1];
    $Y[$i][$j] = $pgap*$Y1 + $pgap*$Y2 + $pgap*$Y3; 
    $X1 = $MX*$M[$i-1][$j]; 
    $X2 = $XX*$X[$i-1][$j];
    $X3 = $BX*$B[$i-1][$j];
    $X[$i][$j] = $pgap*$X1 + $pgap*$X2 + $pgap*$X3;
    
    $f[$i][$j] = $M[$i][$j] + $X[$i][$j] + $Y[$i][$j];
    print("i=$i j=$j M[$i][$j]=$M[$i][$j] X[$i][$j]=$X[$i][$j] Y[$i][$j]=$Y[$i][$j]");
    print"\n";
  }
}

#########################  Backward Probability ##########################################

#### Backward initialization ####
$x = length($seq2); 
$y = length($seq1); 

#### Transition probabilities ####
$BX = $delta; $BY = $delta; $BM = 1 - 2*$delta - $tau; 
#$YB = $delta; $XB = $delta; $MB = 1 - 2*$delta - $tau;
$XM = $delta; $YM = $delta; $MM = 1-2*$delta - $tau;  
$XX = $eta; $XY = 0; $MX = 1 - $eta - $tau; 
$YY = $eta; $YX = 0; $MY = 1 - $eta - $tau;

#### Backward Initialization ####
$B[$x][$y] = 1; $B[$x][$y-1] = 0; $B[$x-1][$y] = 0;
$M[$x][$y] = 0; $M[$x][$y-1] = 0; $M[$x-1][$y] = 0;
$X[$x-1][$y] = $pgap*$BX*$B[$x][$y];
$Y[$x][$y-1] = $pgap*$BY*$B[$x][$y];
$X[$x][$y] = 0; $X[$x][$y-1] = 0;
$Y[$x][$y] = 0; $Y[$x-1][$y] = 0; 
$T[$x][$y-1] = 'L'; $T[$x-1][$y] = 'U'; 

#for(my $i=length($seq2); $i>=0; $i--){
#  $M[$i][$y] = 0; $X[$i][$y] = $pgap*$XX*X[$i+1][$y]; 
#  $Y[$i][0] = 0; $T[$i][0] = 'U';  
#}
#for(my $j=length($seq1); $i>=0; $j--){
#  $M[$x][$j] = 0; $X[$x][$j] = 0; 
#  $Y[$x][$j] = $pgap*$YY*Y[$x][$j+1]; $T[$x][$j] = 'L';  
#}

for(my $i=1; $i<=length($seq2); $i++){
  for(my $j=1; $j<=length($seq1); $j++){
    $p[$i][$j] = $f[$i][$j]*$b[$i-1][$j-1]/$f[length($seq2)][length($seq1)]; 
  }
}

#### Initialize M ####
for(my $i=length($seq2)-2; $i>=0; $i--){
  $M[$i][$y] = 0;
}
for(my $j=length($seq1)-2; $j>=0; $j--){
  $M[$x][$j] = 0;
}

#### Initialize X ####
for(my $i=length($seq2)-2; $i>=0; $i--){
  $X[$i][$y] = $pgap*$XX*$X[$i+1][$y];
}
for(my $j=length($seq1)-2; $j>=0; $j--){
  $X[$x][$j] = 0; 
}

#### Initialize Y #### 
for(my $i=length($seq2)-2; $i>=0; $i--){
  $Y[$i][$y] = 0; 
}
for(my $j=length($seq1)-2; $j>=0; $j--){
  $Y[$x][$j] = $pgap*$YY*$Y[$x][$j+1]; 
}

print"\n";
print"\n"; 

#### Backwards Recurrence ####
for(my $i=length($seq2)-1; $i>=0; $i--){
  for(my $j=length($seq1)-1; $j>=0; $j--){
    $ch1 = substr($seq2, $i, 1); 
    $ch2 = substr($seq1, $j, 1);
    if($ch1 eq $ch2){
      $MorMM = $pm;
    }
    else{
      $MorMM = $pmm; 
    }
    $M1 = $MM*$M[$i+1][$j+1];
    $M2 = $MX*$X[$i+1][$j+1];
    $M3 = $YM*$Y[$i+1][$j+1];
    $M4 = $BM*$B[$i+1][$j+1];  
    $M[$i][$j] = $MorMM*$M1 + $MorMM*$M2 + $MorMM*$M3 + $MorMM*$M4;
    $Y1 = $MY*$M[$i][$j+1];
    $Y2 = $YY*$Y[$i][$j+1]; 
    $Y3 = $YB*$B[$i][$j+1]; 
    $Y[$i][$j] = $pgap*$Y1 + $pgap*$Y2 + $pgap*$Y3; 
    $X1 = $XM*$M[$i+1][$j]; 
    $X2 = $XX*$X[$i+1][$j];
    $X3 = $XB*$B[$i+1][$j];
    $X[$i][$j] = $pgap*$X1 + $pgap*$X2 + $pgap*$X3;
    $b[$i][$j] = $M[$i][$j] + $X[$i][$j] + $Y[$i][$j];
    print("i=$i j=$j M[$i][$j]=$M[$i][$j] X[$i][$j]=$X[$i][$j] Y[$i][$j]=$Y[$i][$j]");
    print"\n";
  }
}

$pxy = $f[$x][$y];
print "pxy is: $pxy";

##### Finding Posterior Probability ####
#for(my $i=1; $i<=length($seq2); $i++){
#  for(my $j=1; $j<=length($seq1); $j++){
#    $p[$i][$j] = $f[$i][$j]*$b[$i][$j]/$f[length($seq2)][length($seq1)];
#  }
#}
for(my $i=0;$i<length($seq2);$i++){
  for(my $j=0;$j<length($seq1);$j++){
    $p[$i][$j] = $f[$i+1][$j+1]*$b[$i][$j]; 
    $p[$i][$j] /= $pxy;
  }
}

print"\n"; print"\n";
print "Matrix f is\n"; 
for(my $i=1; $i<=length($seq2); $i++){
  for(my $j=1; $j<=length($seq1); $j++){
    print "$f[$i][$j] "; 
  }
  print "\n"; 
}
print "\n"; 

print "Matrix b is\n"; 
for(my $i=0; $i<length($seq2); $i++){
  for(my $j=0; $j<length($seq1); $j++){
    print "$b[$i][$j] "; 
  }
  print "\n"; 
}
print "\n";

print "Matrix p is\n"; 
for(my $i=0; $i<=length($seq2); $i++){
  for(my $j=0; $j<=length($seq1); $j++){
    print "$p[$i][$j] "; 
  }
  print "\n"; 
}
print "\n";