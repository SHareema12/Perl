##1st argument, file name that we need to run affine gaps on 
$file = shift;
open(IN,$file);
@data = <IN>;
chomp(@data); 

$seq1 = $data[1];
$seq2 = $data[3];
#
#print "$seq1\n";
#print "$seq2\n";

#$match = 5; $mismatch = -4;

##2nd argument: gap open score
$gap = shift;

##3rd argument: gap extension score
$e = shift; 

##4th argument: blosum62 matrix 
##Read BLOSUM62 matrix from file and use those subsitution costs 
##Read BLOSUM matrix from file into a hashtable 
%blosum = (); 
open(IN, "BLOSUM_62.txt"); 
$aa = <IN>; 
@aa = split(/\s+/, $aa); 
while ($l = <IN>){
  chomp $l; 
  @l = split(/\s+/,$l);
  for(my $i=1; $i<scalar(@l); $i++){
    $key = $l[0] . $aa[$i]; 
    $blosum{"$key"} = $l[$i]; 
  }
}

@V = ();
@E = ();
@F = ();
@M = ();
@T = (); 

##Initialization of all matrices. Infinity values are there b/c we can't align a sequence to all gaps
$V[0][0] = 0; $M[0][0] = 0; 
$E[0][0] = -inf; $F[0][0] = -inf; 
$T[0][0] = 0;
for(my $i=1; $i<=length($seq1);$i++){
   $V[$i][0] = $gap + ($i-1)*$e;
   $F[$i][0] = $gap + ($i-1)*$e;
   $M[$i][0] = -inf; 
   $E[$i][0] = -inf;
   $T[$i][0] = "U";  
}
for(my $j=1; $j<=length($seq2);$j++){
   $V[0][$j] = $gap + ($j-1)*$e;
   $E[0][$j] = $gap + ($j-1)*$e;
   $M[0][$j] = -inf;
   $F[0][$j] = -inf; 
   $T[0][$j] = "L"; 
}
#print "M matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($M[$i][$j] , " ");
#  }
#  print("\n");
#}
#
#print "F matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($F[$i][$j] , " ");
#  }
#  print("\n");
#}
#print "E matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($E[$i][$j] , " ");
#  }
#  print("\n");
#}


##Recurrence 
##Fill M matrix first-diagonal of V value+m or mm. Diagonal values of V are already there
###M[i,j) = V{i-1,j-1) + (m+mm)  
for(my $i=1;$i<=length($seq1);$i++){
   for(my $j=1;$j<=length($seq2);$j++){
      $ch1 = substr($seq1,$i-1,1);
      $ch2 = substr($seq2,$j-1,1);
#      if($ch1 eq $ch2){
#          $MorMM =  $match; 
#      }
#      else{
#          $MorMM = $mismatch;
#      }
      $MorMM = $blosum{substr($seq1, $i-1, 1).substr($seq2,$j-1,1)};
      ##Filling in M matrix 
      $M[$i][$j] = $V[$i-1][$j-1] + $MorMM; 
      $diag = $M[$i][$j]; 
      ##Filling in F matrix
      $MvalueF = $M[$i-1][$j] + $gap;
      $MvalueE = $M[$i][$j-1] + $gap; 
      $Evalue = $E[$i][$j-1] + $e;
      $Fvalue = $F[$i-1][$j] + $e; 
      if($MvalueF >= $Fvalue){
         $F[$i][$j] = $MvalueF;
      }  
      else{
         $F[$i][$j] = $Fvalue;
      }
      $up = $F[$i][$j];
      if($MvalueE >= $Evalue){
         $E[$i][$j] = $MvalueE; 
      }  
      else{
         $E[$i][$j] = $Evalue;
      }
      $left = $E[$i][$j];
      ##Fill in T and V matrix. 'D'-->if V=M, 'L' --> if V=E, 'U' --> if V=F
      ## V[i][j] = max(E[i][j],F[i][j].M[i][j])
      if($up >= $diag && $up >= $left){
      $V[$i][$j] = $up;
      $T[$i][$j] = 'U';
      }
      elsif($left >= $diag && $left >= $up){
      $V[$i][$j] = $left;
      $T[$i][$j] = 'L';
      }
      else{
      $V[$i][$j] = $diag;
      $T[$i][$j] = 'D';
      }
   } 
}
#print "M matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($M[$i][$j] , " ");
#  }
#  print("\n");
#}
#
#print "F matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($F[$i][$j] , " ");
#  }
#  print("\n");
#}
#print "E matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($E[$i][$j] , " ");
#  }
#  print("\n");
#}
#print "V matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($V[$i][$j] , " ");
#  }
#  print("\n");
#}
#print "T matrix\n"; 
#for(my $i=0; $i<=length($seq1);$i++){
#  for(my $j=0; $j<=length($seq2);$j++){
#    print($T[$i][$j] , " ");
#  }
#  print("\n");
#}


####Traceback 
$aligned_seq1 = "";
$aligned_seq2 = "";
$i = length($seq1);
$j = length($seq2);
while(!($i == 0 && $j == 0)){
	if($T[$i][$j] eq "L"){
    $ch1 = '-'; 
    $ch2 = substr($seq2, $j-1, 1); 
    $j--;
   }
  elsif($T[$i][$j] eq "U"){
    $ch1 = substr($seq1,$i-1,1);
    $ch2 = '-';
    $i--;
  } 
  else{
    $ch1 = substr($seq1,$i-1,1);
    $ch2 = substr($seq2,$j-1,1);
  $i--;
  $j--;
  }
  $aligned_seq1 = $ch1 . $aligned_seq1; 
  $aligned_seq2 = $ch2 . $aligned_seq2; 
}

#print "$data[0]\n"; 
#print "$aligned_seq1\n";
#print "$data[2]\n"; 
#print "$aligned_seq2\n"; 
#
#print "$V[length($seq1)][length($seq2)]/n";

#//
#//print "$data[0]\n"; 
#//print "$aligned_seq1\n";
#//print "$data[2]\n"; 
#//print "$aligned_seq2\n"; 

$tempfile = "AiNWAGBL.txt";

if(-f $tempfile){
      unlink $tempfile;
    }

open(OUTFILE,">>$tempfile") or die;
print OUTFILE "$data[0]\n"; 
print OUTFILE "$aligned_seq1\n";
print OUTFILE "$data[2]\n"; 
print OUTFILE "$aligned_seq2\n"; 
close(OUTFILE); 
