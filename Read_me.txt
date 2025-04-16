This is the last version of script of TV-PRO-seq data processing:

1.If you applied single end sequencing , the data could be directly used for analysis, if paired end sequencing was used, only use R1 files for further analysis.
2.Move all the .fastq/.fastq.gz files into a new folder, download the PRO_raw_processing.sh script from github to the same folder. Change the BOWTIE2_INDEX with your path then directly run ./PRO_raw_processing.sh which process .fastq/.fastq.gz file to bedgraph file.
3.Merge the bedgraph files into two files, arrange the bedgraph file as follow:
perl bedgraph_merge.pl 0.5_1_p.bedgraph 0.5_2_p.bedgraph 2_1_p.bedgraph 2_2_p.bedgraph 8_1_p.bedgraph 8_2_p.bedgraph 32_1_p.bedgraph 32_2_p.bedgraph > merge_p.bedgraph
perl bedgraph_merge.pl 0.5_1_m.bedgraph 0.5_2_m.bedgraph 2_1_m.bedgraph 2_2_m.bedgraph 8_1_m.bedgraph 8_2_m.bedgraph 32_1_m.bedgraph 32_2_m.bedgraph > merge_m.bedgraph
4.Preform peak calling as follow:
perl TV_callpeak.pl merge_p.bedGraph merge_m.bedGraph
This script will generate two files, TVPRO_peak file have record original reads of all peaks13, chr_read file record the total reads of mitochondrion and nucleus.
5.Generate the final result by command line: python average_time_count.py; The final pausing time estimation is inside file pausing_times.
