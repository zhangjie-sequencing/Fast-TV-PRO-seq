#!/usr/bin/perl
# This script is use for extract peaks from Merged bedgraph file from 4 time point of TVPRO-seq, input should be merged data of plus strand then minus strand

# perl TV_callpeak.pl merge_p.bedGraph merge_m.bedGraph

# create output file
open ( my $outfile, ">TVPRO_peak" ) or die; #data output
print $outfile "chr strand position total_reads read_0.5m read_2m read_8m read_32m surround7 surround201\n";
my @empty_array;
for (my $var = 0; $var < 201; $var++) {
	push @empty_array,0;
}

# reading input file

for (my $file = 0; $file < 2; $file++) {
	my @Temp_lines=@empty_array;
	my @Temp_reads=@empty_array;
	my @Temp_pos=@empty_array;
	my @Min1=@empty_array;

	open FILE, $ARGV[$file];
	my $title=<FILE>;
	my $title2=<FILE>;
	my $strand;
	my $inputname=$ARGV[$file];
	if ($file eq 0) {
		$strand="+";
	}else{
		$strand="-";
	}
	print "Now analysis $strand strand file $inputname\n";
	# Start peak selection:
	# peak should be the highest one of +-3 region and bigger than 12, also if 
	my $last_chr="none";
	my $counter=1;
	while (defined (my $line = <FILE>) ) {
		chomp $line;
		my @parameters=split(/ /,$line);
		my $chr=$parameters[0];
		my $pos=$parameters[1];
		my $reads=$parameters[2];
		$counter++;
		if ($counter % 1000000 eq 0) {
			print "Now analysising $counter line, $chr $pos\n";
		}
		my $min1_check;
		
		$T05=$parameters[3]+$parameters[4];
        $T2=$parameters[5]+$parameters[6];
        $T8=$parameters[7]+$parameters[8];
        $T32=$parameters[9]+$parameters[10];

		if ($T05>0 and $T2>0 and $T8>0 and $T32>0) {
			$min1_check=1;
		}else{
			$min1_check=0;
		}
		if ($chr eq $last_chr) {
			$chr_counter++;
			
			# Loading peak in to template array and analysis the line 100 before current line
			push @Min1,$min1_check;
			push @Temp_lines,"$T05 $T2 $T8 $T32";
			push @Temp_reads,$reads;
			push @Temp_pos,$pos;
			shift @Min1;
			shift @Temp_lines;
			shift @Temp_reads;
			shift @Temp_pos;
			if ($Temp_reads[100]>12) {
				my $check=0;
				for (my $i = 97; $i < 104; $i++) {
					if (abs($Temp_pos[100]-$Temp_pos[$i])<4){
						if ($Temp_reads[$i]>$Temp_reads[100]) {
							$check++;
						}
					}	
				}
				if ($check eq 0) {
					my $surround7=0;
					my $surround201=0;
					for (my $i = 0; $i < 201; $i++) {
						if (abs($Temp_pos[100]-$Temp_pos[$i])<101) {
							$surround201=$surround201+$Temp_reads[$i];
						}
					}
					for (my $i = 97; $i < 104; $i++) {
						if (abs($Temp_pos[100]-$Temp_pos[$i])<4) {
							$surround7=$surround7+$Temp_reads[$i];
						}
					}
					if ($Min1[100] eq 1) {
						print $outfile "$chr $strand $Temp_pos[100] $Temp_reads[100] $Temp_lines[100] $surround7 $surround201\n";
					}
				}	
			}
		}else{
			# Output the peaks have not been count and reset the arrays
			for (my $scan = 101; $scan < 201; $scan++) {
				if ($Temp_reads[$scan]>12) {
					my $check;
					for (my $var = 0; $var < 201; $var++) {
						if ($Temp_reads[$var]>$Temp_reads[$scan] and abs($Temp_pos[$var]-$Temp_pos[$scan]<4)) {
							$check++;
						}
					}
					if ($check eq 0) {
						my $surround7=0;
						my $surround201=0;
						for (my $var = 0; $var < 201; $var++) {
							if (abs($Temp_pos[$var]-$Temp_pos[$scan]<4)) {
								$surround7=$surround7+$Temp_pos[$var];
							}
							if (abs($Temp_pos[$var]-$Temp_pos[$scan]<101)) {
								$surround201=$surround201+$Temp_pos[$var];
							}
						}
						if ($Min1[$scan] eq 1) {
						print $outfile "$chr $strand $Temp_pos[$scan] $Temp_reads[$scan] $Temp_lines[$scan] $surround7 $surround201\n";
						}
					}
				}
			}
			# reset virables:
			
			$last_chr=$chr;
			@Min1=@empty_array;
			@Temp_lines=@empty_array;
			@Temp_reads=@empty_array;
			@Temp_pos=@empty_array;
			push @Min1,$min1_check;
			push @Temp_lines,"$T05 $T2 $T8 $T32";
			push @Temp_reads,$reads;
			push @Temp_pos,$pos;
			shift @Min1;
			shift @Temp_lines;
			shift @Temp_reads;
			shift @Temp_pos;
		}
	}
}






print "Now generating chr_read file\n";
my @sum_m = (0) x 8;      
my @sum_other = (0) x 8;  

foreach my $file (@ARGV) {  
    open(my $fh, '<', $file) or die "could not open $file: $!";
    <$fh>; <$fh>;  
    
    while (my $line = <$fh>) {
        chomp $line;
        my @cols = split(/\s+/, $line);
        my $chr = $cols[0];
        

        if ($chr =~ /^(M|chrM|MtDNA)$/) {
            for my $i (3..10) {  
                $sum_m[$i-3] += ($cols[$i] || 0);
            }
        } else {
            for my $i (3..10) {
                $sum_other[$i-3] += ($cols[$i] || 0);
            }
        }
    }
    close $fh;
}


sub merge_columns {
    my @arr = @_;
    my @merged;
    for (my $i = 0; $i < @arr; $i += 2) {
        push @merged, ($arr[$i] + $arr[$i+1]);
    }
    return @merged;
}

my @merged_m = merge_columns(@sum_m);
my @merged_other = merge_columns(@sum_other);

open(my $out, '>', 'chr_read') or die "无法创建chr_read文件: $!";
print $out join(" ", @merged_m) . "\n";
print $out join(" ", @merged_other) . "\n";
close $out;