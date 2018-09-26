### All SAM files have been generated and are in the directory ###
### This program determines accuracy of each of the different reads ### 

### Take input of SAM file ### 
$file = shift; 
$read = $file;
print"Read: $read\n";
open(IN, $file);
@data = <IN>;
chomp(@data);

$columned_data = (); 

## Split data and put into 2d array, only taking first 4 columns ##
for(my $i=2;$i<@data;$i++){
  @array1 = split /\s+/, $data[$i];
  $x = @array1; 
#  print"column size for row $i: $x\n";
  for(my $j =0;$j<5;$j++){
    $columned_data[$i-2][$j] = $array1[$j];
  } 
}

$total_alignments = @columned_data; 
#for my $row (@columned_data) {
#    print join(",", @{$row}), "\n";
#}

### Taking substring of positions to compare alignment ###
for (my $i=0;$i<@columned_data;$i++){
  $columned_data[$i][0] = substr($columned_data[$i][0],6);
  $columned_data[$i][0] = substr($columned_data[$i][0],0,(length($columned_data[$i][0]))-2);  
}

#for my $row (@columned_data) {
#    print join(",", @{$row}), "\n";
#}

### Comparing alignment ###
$aligned_correct = 0; 
$aligned_count = 0;
for (my $i=0;$i<@columned_data;$i++){
  if($columned_data[$i][1] eq "16" || $columned_data[$i][1] eq "0"){
    $aligned_count += 1; 
    $true_pos = int($columned_data[$i][0]); 
    $aligned_pos = int($columned_data[$i][3]); 
    #print "Alignment: $i true position: $true_pos aligned position: $aligned_pos\n";
    $diff = $aligned_pos - $true_pos; 
    if($diff >= -100 && $diff <= 100){
      $aligned_correct += 1; 
    }
  }
}

$accuracy = ($aligned_correct/$total_alignments)*100;
print"\n";
print"Read: $read\n";
print"Total number aligned: $aligned_count\n"; 
print"Total number correctly aligned: $aligned_correct\n";
print"Total number of reads: $total_alignments\n"; 
print"Accuracy: $accuracy%\n";
print"\n";


 