$true_alignments_file = shift; 
open(TRUE, $true_alignments_file);
@traindata = <TRUE>;
chomp(@traindata);

##Gap open values to test
@go = (-16, -14, -12, -10, -8, -6, -4, -2);

##Gap extension values to test
@ge = (-2, -1, -.5, -.1);

%accuracy = (); 
$numAlignments = scalar(@traindata); 
print "Num alignments: $numAlignments\n";

for(my $i=0;$i<@traindata;$i++){
  $unaligned = $traindata[$i];
  $aligned = $unaligned; 
  $aligned =~ s/\.unaligned//;
  for (my $j=0;$j<@go;$j++){ 
    $gopen = $go[$j]; 
     for (my $k=0;$k<@ge;$k++){
        $gext = $ge[$k];
        $hashkey = $gopen . '/' . $gext;   

        ##Running affine gaps on unaligned sequence with test gap open and gap extension 
        $agCall = ('perl affinegaps.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $agCall .= $unaligned . '" ' . $gopen . " " . $gext;
        system($cmd); 
        
        ##Result from test affine gaps is saved in temp file AiAG, so now run accuracy b/w ALIGNED and            ##alignment we got with our test gap open and gap extension
        $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $accCall .= $aligned . '" AiNWAG.txt'; 
        system($cmd);
        
        ##Retrieve accuracy from above system call 
        open(IN,"acc.txt");
        @acc = <IN>; 
        chomp(@acc);
        close(IN);
        $tempAcc = $acc[0]; 
        print "Gap: $gopen and ext: $gext Accuracy: $acc[0]\n";  
        
        ##Add the accurac to the hash table with corresponding gap open and gap extension value 
        $accuracy{$hashkey} = $accuracy{$hashkey} + $tempAcc;          
     }
  }
}

##Get keys of accuracy table 
@acc_keys = keys(%accuracy); 
$maxAcc = 0;
$bestGopen = 0; 
$bestGext = 0;  
$len = scalar(@acc_keys);
for (my $i=0; $i<scalar(@go); $i++){
  for(my $j=0; $j<scalar(@ge); $j++){
    $accKey = $go[$i] . '/' . $ge[$j]; 
    $b4avg = $accuracy{$accKey}; 
    #print "Before averaging: $b4avg\n";
    $afteravg = $b4avg/$numAlignments;
    $accuracy{$accKey} = $afteravg; 
    $afteravg = $accuracy{$accKey}; 
    #print "After averaging: $afteravg\n";
    if ($accuracy{$accKey} >= $maxAcc){
      $maxAcc = $accuracy{$accKey}; 
      $bestGopen = $go[$i]; 
      $bestGext = $ge[$j]; 
    }
  }  
}
$lowestError = 1-$maxAcc; 
print "Best Gap Open: $bestGopen\n"; 
print "Best Gap Extension: $bestGext\n"; 
print "Best Accuracy:$maxAcc\n"; 
print "Lowest Error: $lowestError\n"; 


