#!/bin/bash

# Initialize variables
vcf=""
RegeniefileName=""
outputDir=""

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --vcf) vcf="$2"; shift ;;  # Get the value for --vcf
        --RegeniefileName) RegeniefileName="$2"; shift ;; # Get the value for --RegeniefileName
        --outputDir) outputDir="$2"; shift ;; # Get the value for --outputDir
        *) echo "Unknown parameter passed: $1"; exit 1 ;; # Catch unknown parameters
    esac
    shift
done

# Check if the required arguments are provided
if [[ -z "$vcf" || -z "$RegeniefileName" || -z "$outputDir" ]]; then
    echo "Error: --vcf, --RegeniefileName, and --outputDir are required."
    exit 1
fi


# Submit jobs with different annotations, passing vcf and filename as arguments
sbatch step9.sh "list_annot_run1" Run1_Count  "$vcf" "$RegeniefileName" "$outputDir"
sbatch step9.sh "list_annot_run2" Run2_Count  "$vcf" "$RegeniefileName" "$outputDir"
sbatch step9.sh "list_annot_run3" Run3_Count  "$vcf" "$RegeniefileName" "$outputDir"
sbatch step9.sh "list_annot_run4" Run4_Count  "$vcf" "$RegeniefileName" "$outputDir"
