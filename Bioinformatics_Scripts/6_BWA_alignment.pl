### This script runs the bwa_alignment_script for each of the 6 reads 

$cmd = 'perl bwa_alignment_script.pl reads_01.sam > BWAres01'; 
system($cmd); 

$cmd = 'perl bwa_alignment_script.pl reads_05.sam > BWAres05'; 
system($cmd);

$cmd = 'perl bwa_alignment_script.pl reads_1.sam > BWAres1'; 
system($cmd);

$cmd = 'perl bwa_alignment_script.pl reads_15.sam > BWAres15'; 
system($cmd);

$cmd = 'perl bwa_alignment_script.pl reads_2.sam > BWAres2'; 
system($cmd);

$cmd = 'perl bwa_alignment_script.pl reads_25.sam > BWAres25'; 
system($cmd);

#### Open up the result and print out alignment scores for each read ####
open(IN, "BWAres01");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}

open(IN, "BWAres05");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}

open(IN, "BWAres1");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}

open(IN, "BWAres15");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}

open(IN, "BWAres2");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}

open(IN, "BWAres25");
@data = <IN>;
chomp(@data);

for (my $i=0;$i<@data;$i++){
  print"$data[$i]\n";
}