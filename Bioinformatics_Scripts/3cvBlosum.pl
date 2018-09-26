$true_alignments_file = shift; 
open(TRUE, $true_alignments_file);
@traindata = <TRUE>;
chomp(@traindata);

##Gap open values to test
@go = (-16, -14, -12, -10, -8, -6, -4, -2);

##Gap extension values to test
@ge = (-2, -1, -.5, -.1);

%accuracyNW = ();
%accuracyNWBL = ();
%accuracyNWAG = (); 
%accuracyNWAGBL = ();  
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

        ##Running normal affine gaps on unaligned sequence with test gap open and gap extension 
        $agCall = ('perl affinegaps.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $agCall .= $unaligned . '" ' . $gopen . " " . $gext;
        print "system command: $cmd\n"; 
        system($cmd); 
        
        ##Result from test affine gaps is saved in temp file AiAG, so now run accuracy b/w ALIGNED and            ##alignment we           got with our test gap open and gap extension
        $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $accCall .= $aligned . '" AiNWAG.txt'; 
        system($cmd);
        
        ##Retrieve accuracy from above system call 
        open(IN,"acc.txt");
        @acc = <IN>; 
        chomp(@acc);
        close(IN);
        $tempAcc = $acc[0]; 
        print "AG Gap: $gopen and ext: $gext Accuracy: $acc[0]\n";  
        
        ##Add the accurac to the hash table with corresponding gap open and gap extension value 
        $accuracyNWAG{$hashkey} = $accuracyNWAG{$hashkey} + $tempAcc;     
        
        ##Running affine gaps with BLOSUM 62 scores on unaligned sequence with test gap open and gap          extension 
        $agCall = ('perl affinegaps_blosum62.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $agCall .= $unaligned . '" ' . $gopen . " " . $gext;
        system($cmd); 
        
        ##Result from test affine gaps is saved in temp file AiAGBL, so now run accuracy b/w ALIGNED          and ##alignment we got with our test gap open and gap extension
        $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
        $cmd = $accCall .= $aligned . '" AiNWAGBL.txt'; 
        system($cmd);
        
        ##Retrieve accuracy from above system call 
        open(IN,"acc.txt");
        @acc = <IN>; 
        chomp(@acc);
        close(IN);
        $tempAcc = $acc[0]; 
        print "AGBL Gap: $gopen and ext: $gext Accuracy: $acc[0]\n";  
        
        ##Add the accurac to the hash table with corresponding gap open 
        $accuracyNWAGBL{$hashkey} = $accuracyNWAGBL{$hashkey} + $tempAcc;         
     }
    $hashkey = $gopen;

    ##Running normal NW on unaligned sequence with test gap open and gap extension 
    $nwCall = ('perl nw.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
    $cmd = $nwCall .= $unaligned . '" ' . $gopen;
    system($cmd); 
        
    ##Result from test NW is saved in temp file AiNW, so now run accuracy b/w ALIGNED and            #    #alignment we got with our test gap open and gap extension
    $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
    $cmd = $accCall .= $aligned . '" AiNW.txt'; 
    system($cmd);
        
    ##Retrieve accuracy from above system call 
    open(IN,"acc.txt");
    @acc = <IN>; 
    chomp(@acc);
    close(IN);
    $tempAcc = $acc[0]; 
    print "NWGap: $gopen Accuracy: $acc[0]\n";  
       
    ##Add the accuracy to the hash table with corresponding gap open value 
    $accuracyNW{$hashkey} = $accuracyNW{$hashkey} + $tempAcc; 
        
    ##Running NW with BLOSUM 62 scores on unaligned sequence with test gap open 
    $nwCall = ('perl nw_blosum62.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
    $cmd = $nwCall .= $unaligned . '" ' . $gopen;
    system($cmd); 
        
    ##Result from test affine gaps is saved in temp file AiNWBL, so now run accuracy b/w ALIGNED and       ##alignment we got with our test gap open and gap extension
    $accCall = ('perl accuracy.pl "/afs/cad.njit.edu/courses/bnfo/f17/bnfo/601/001/sh527/');
    $cmd = $accCall .= $aligned . '" AiNWBL.txt'; 
    system($cmd);
        
    ##Retrieve accuracy from above system call 
    open(IN,"acc.txt");
    @acc = <IN>; 
    chomp(@acc);
    close(IN);
    $tempAcc = $acc[0]; 
    print "NWBL Gap: $gopen Accuracy: $acc[0]\n";  
       
    ##Add the accuracy to the hash table with corresponding gap open and gap extension value 
    $accuracyNWBL{$hashkey} = $accuracyNWBL{$hashkey} + $tempAcc;         
  }
}

##Get keys of accuracy table 
@accNW_keys = keys(%accuracyNW); 
@accNWBL_keys = keys(%accuracyNWBL);
@accNWAG_keys = keys(%accuracyNWAG);
@accNWAGBL_keys = keys(%accuracyNWAGBL);
$maxAccNWAG = 0;
$bestGopenNWAG = 0; 
$bestGextNWAG = 0; 
$maxAccNWAGBL = 0;
$bestGopenNWAGBL = 0; 
$bestGextNWAGBL = 0; 
$maxAccNW = 0;
$bestGopenNW = 0; 
$maxAccNWBL = 0;
$bestGopenNWBL = 0; 
$len = scalar(@accNWAG_keys);
for (my $i=0; $i<scalar(@go); $i++){
  for(my $j=0; $j<scalar(@ge); $j++){
    $accKey = $go[$i] . '/' . $ge[$j]; 
    $b4avg = $accuracyNWAG{$accKey}; 
    #print "Before averaging: $b4avg\n";
    $afteravg = $b4avg/$numAlignments;
    $accuracyNWAG{$accKey} = $afteravg; 
    $afteravg = $accuracyNWAG{$accKey}; 
    #print "After averaging: $afteravg\n";
    if ($accuracyNWAG{$accKey} >= $maxAccNWAG){
      $maxAccNWAG = $accuracyNWAG{$accKey}; 
      $bestGopenNWAG = $go[$i]; 
      $bestGextNWAG = $ge[$j]; 
    }
    $b4avg = $accuracyNWAGBL{$accKey}; 
    #print "Before averaging: $b4avg\n";
    $afteravg = $b4avg/$numAlignments;
    $accuracyNWAGBL{$accKey} = $afteravg; 
    $afteravg = $accuracyNWAGBL{$accKey}; 
    #print "After averaging: $afteravg\n";
    if ($accuracyNWAGBL{$accKey} >= $maxAccNWAGBL){
      $maxAccNWAGBL = $accuracyNWAGBL{$accKey}; 
      $bestGopenNWAGBL = $go[$i]; 
      $bestGextNWAGBL = $ge[$j]; 
    }
  } 
  $accKey = $go[$i]; 
  $b4avg = $accuracyNW{$accKey}; 
  #print "Before averaging: $b4avg\n";
  $afteravg = $b4avg/$numAlignments;
  $accuracyNW{$accKey} = $afteravg; 
  $afteravg = $accuracyNW{$accKey}; 
  #print "After averaging: $afteravg\n";
  if ($accuracyNW{$accKey} >= $maxAccNW){
    $maxAccNW = $accuracyNW{$accKey}; 
    $bestGopenNW = $go[$i]; 
  } 
  
  $b4avg = $accuracyNWBL{$accKey}; 
  #print "Before averaging: $b4avg\n";
  $afteravg = $b4avg/$numAlignments;
  $accuracyNWBL{$accKey} = $afteravg; 
  $afteravg = $accuracyNWBL{$accKey}; 
  #print "After averaging: $afteravg\n";
  if ($accuracyNWBL{$accKey} >= $maxAccNWBL){
    $maxAccNWBL = $accuracyNWBL{$accKey}; 
    $bestGopenNWBL = $go[$i]; 
  }
}
$lowestErrorNWAG = 1-$maxAccNWAG; 
$lowestErrorNWAGBL = 1-$maxAccNWAGBL; 
$lowestErrorNW = 1-$maxAccNW; 
$lowestErrorNWBL = 1-$maxAccNWBL; 

print "Needleman-Wunsch\n";
print "Best Gap Open: $bestGopenNW\n";  
print "Best Accuracy:$maxAccNW\n"; 
print "Lowest Error: $lowestErrorNW\n"; 

print "Needleman-Wunsch BLOSUM62\n";
print "Best Gap Open: $bestGopenNWBL\n";  
print "Best Accuracy:$maxAccNWBL\n"; 
print "Lowest Error: $lowestErrorNWBL\n"; 

print "Needleman-Wunsch Affine Gaps\n";
print "Best Gap Open: $bestGopenNWAG\n"; 
print "Best Gap Extension: $bestGextNWAG\n"; 
print "Best Accuracy:$maxAccNWAG\n"; 
print "Lowest Error: $lowestErrorNWAG\n"; 

print "Needleman-Wunsch Affine Gaps BLOSUM62\n";
print "Best Gap Open: $bestGopenNWAGBL\n"; 
print "Best Gap Extension: $bestGextNWAGBL\n"; 
print "Best Accuracy:$maxAccNWAGBL\n"; 
print "Lowest Error: $lowestErrorNWAGBL\n"; 