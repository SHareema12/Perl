### Input: X and Y, going to hash both based on kmers and increment if they match ###
$file = shift; 
open(IN, $file);
@X_data = <IN>; 
chomp(@X_data);

$file1 = shift; 
open(IN, $file1);
@Y_data = <IN>;
chomp(@Y_data); 

$seed = shift; ##seed I guess 

$X = $X_data[0];
$Y = $Y_data[0]; 

%hash_X = (); 
$count = (length($X)-length($seed))+1; 
$count2 = (length($Y)-length($seed))+1; 
for (my $i = 0; $i<$count; $i++){
  $key = ""; 
  for(my $j=0; $j<length($seed);$j++){
    
  }
}