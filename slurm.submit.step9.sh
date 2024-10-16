#!/bin/bash

sbatch step9.sh "list_annot_run1" Run1_Count  $1 $2
sbatch step9.sh "list_annot_run2" Run2_Count  $1 $2
sbatch step9.sh "list_annot_run3" Run3_Count  $1 $2
sbatch step9.sh "list_annot_run4" Run4_Count  $1 $2
