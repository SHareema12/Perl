$true_alignments_file = shift; 
open(TRUE, $true_alignments_file);
@traindata = <TRUE>;
chomp(@traindata);

#### delta and eta values to test ####
@delta = (.1,.15,.2,.25,.3,.35,.4,.45,.5); 
@eta = (.2,.25,.3,.35,.4,.45,.5,.55,.6);

%accuracy = ();
$numAlignments = scalar(@traindata);
for(my $i=0;$i<@traindata;$i++){
  $unaligned = $traindata[$i];
  $aligned = $unaligned; 
  $aligned =~ s/\.unaligned//;
  for (my $j=0;$j<@delta;$j++){ 
    $deltaval = $delta[$j]; 
     for (my $k=0;$k<@eta;$k++){
       $etaval = $eta[$k]; 
       $hashkey = $deltaval . '/' . $etaval;
       
        #### Running HMM with current delta and eta values #### 
        $HMMCall = ('perl HMM.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $HMMCall .= $unaligned . '" ' . $deltaval . " " . $etaval;
        print "system command: $cmd\n"; 
        system($cmd); 
        
        ##Result from HMM is saved in file HMMres, so now run accuracy b/w HMM and true sequence
        $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $accCall .= $aligned . '" HMMres.txt'; 
        system($cmd);
        
        ##Retrieve accuracy from above system call 
        open(IN,"acc.txt");
        @acc = <IN>; 
        chomp(@acc);
        close(IN);
        $tempAcc = $acc[0]; 
        print "Delta $delta and eta $eta accuracy: $acc[0]\n";  
        
        ##Add the accurac to the hash table with corresponding delta value 
        $accuracy{$hashkey} = $accuracy{$hashkey} + $tempAcc; 
     }
  }
}

##Get keys of accuracy table 
@acc_keys = keys(%accuracy);
$bestDelta = 0;
$bestEta = 0; 
$maxAcc = 0; 

for (my $i=0; $i<scalar(@delta); $i++){
  for(my $j=0; $j<scalar(@eta); $j++){
    $accKey = $delta[$i] . '/' . $eta[$j]; 
    $b4avg = $accuracy{$accKey}; 
    $accuracy{$accKey} = $accuracy{$accKey}/$numAlignments;; 
    $afteravg = $accuracy{$accKey}; 
    #print "After averaging: $afteravg\n";
    if ($accuracy{$accKey} >= $maxAcc){
      $maxAcc = $accuracy{$accKey}; 
      $bestDelta = $delta[$i]; 
      $bestEta = $eta[$j]; 
    }
  }  
}

foreach my $key (sort { $accuracy{$b} <=> $accuracy{$a} } keys %accuracy) {
    printf "%-8s %s\n", $key, $accuracy{$key};
}

$lowestError = 1-$maxAcc; 
print "Best Delta: $bestDelta\n"; 
print "Best Eta: $bestEta\n"; 
print "Best Accuracy:$maxAcc\n"; 
print "Lowest Error: $lowestError\n"; 


