#!/bin/bash

mkdir -p 01_trimmed_reads 02_fastqc_reports

SUMMARY_CSV="02_fastqc_reports/fastp_summary.csv"
echo -e "Sample\tTotal_Reads\tPassed_Reads\tQ20_Rate\tQ30_Rate\tGC_Content" > "$SUMMARY_CSV"

for R1 in 00_raw_reads/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
  R2="00_raw_reads/${SAMPLE}_R2_001.fastq.gz"

  if [[ -f "$R2" ]]; then
    echo "Processing $SAMPLE..."

    fastp -i "$R1" -I "$R2" \
      -o "01_trimmed_reads/${SAMPLE}_R1.trimmed.fastq.gz" \
      -O "01_trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz" \
      -q 20 -u 30 -n 5 -w 4 \
      -h "02_fastqc_reports/${SAMPLE}_fastp.html" \
      -j "02_fastqc_reports/${SAMPLE}_fastp.json"

    total=$(jq '.summary.before_filtering.total_reads' "02_fastqc_reports/${SAMPLE}_fastp.json")
    passed=$(jq '.summary.after_filtering.total_reads' "02_fastqc_reports/${SAMPLE}_fastp.json")
    q20=$(jq '.summary.after_filtering.q20_rate' "02_fastqc_reports/${SAMPLE}_fastp.json")
    q30=$(jq '.summary.after_filtering.q30_rate' "02_fastqc_reports/${SAMPLE}_fastp.json")
    gc=$(jq '.summary.after_filtering.gc_content' "02_fastqc_reports/${SAMPLE}_fastp.json")

    echo -e "$SAMPLE\t$total\t$passed\t$q20\t$q30\t$gc" >> "$SUMMARY_CSV"
  fi
done
