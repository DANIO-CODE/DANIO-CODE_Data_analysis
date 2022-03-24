#!/bin/bash

sbatch -p high -x kraken --export=ALL -o dome_all_v3.out -J super --mem=8G  --wrap "/usr/bin/python ROSE_main.py -g danrer10 -i ~/projects/DANIO-CODE_Data_analysis/data/SR136_k27acpeaks/dome.gff -r ~/projects/DANIO-CODE_Data_analysis/data/h3k27ac_tagAlign_ln/dome.bam -o dome_all_v3 -t 500"
sbatch -p high -x kraken --export=ALL -o Epi75_all_v3.out -J super --mem=8G --wrap "/usr/bin/python ROSE_main.py -g danrer10 -i ~/projects/DANIO-CODE_Data_analysis/data/SR136_k27acpeaks/Epi75.gff -r ~/projects/DANIO-CODE_Data_analysis/data/h3k27ac_tagAlign_ln/Epi75.bam -o Epi75_all_v3 -t 500"
sbatch -p high -x kraken --export=ALL -o Som8_all_v3.out -J super --mem=8G --wrap "/usr/bin/python ROSE_main.py -g danrer10 -i ~/projects/DANIO-CODE_Data_analysis/data/SR136_k27acpeaks/Som8.gff -r ~/projects/DANIO-CODE_Data_analysis/data/h3k27ac_tagAlign_ln/Som8.bam -o Som8_all_v3 -t 500"
sbatch -p high -x kraken --export=ALL -o Prim5_all_v3.out -J super --mem=8G --wrap "/usr/bin/python ROSE_main.py -g danrer10 -i ~/projects/DANIO-CODE_Data_analysis/data/SR136_k27acpeaks/Prim5.gff -r ~/projects/DANIO-CODE_Data_analysis/data/h3k27ac_tagAlign_ln/Prim5.bam -o Prim5_all_v3 -t 500"
sbatch -p high -x kraken --export=ALL -o LongPec_all_v3.out -J super --mem=8G --wrap "/usr/bin/python ROSE_main.py -g danrer10 -i ~/projects/DANIO-CODE_Data_analysis/data/SR136_k27acpeaks/LongPec.gff -r ~/projects/DANIO-CODE_Data_analysis/data/h3k27ac_tagAlign_ln/LongPec.bam -o LongPec_all_v3 -t 500"

