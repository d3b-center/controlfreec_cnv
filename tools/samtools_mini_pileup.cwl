cwlVersion: v1.0
class: CommandLineTool
id: samtools_mini_pileup
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bvcftools:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: $(inputs.threads)
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      cut -f 1 $(inputs.subset_fai.path) > chr_names.txt

      cut -f 1,2 $(inputs.subset_fai.path) > bed.genome
      
      zcat $(inputs.snp_vcf.path) |  grep -v "#" | awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed

      bedtools sort -i snps.bed -faidx chr_names.txt > $(inputs.snp_vcf.nameroot).safety_sorted.bed

      cat chr_names.txt 
      | xargs -IXX -P $(inputs.threads) sh -c  "samtools mpileup -d 250 -f $(inputs.reference.path) -r XX $(inputs.input_reads.path)
      | awk '{if(\$4 != 0) print \$0}
      | awk {'printf (\"%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", \$1,\$2-1,\$2,\$3,\$4,\$5,\$6)'}
      | bedtools intersect -a stdin -b $(inputs.snp_vcf.nameroot).safety_sorted.bed -g bed.genome -sorted -wa
      | cut -f 1,3- > XX.pileup"

      cat chr_names.txt | xargs -IXX cat XX.pileup >> $(inputs.input_reads.basename).mini.pileup

inputs:
  input_reads: {type: File, secondaryFiles: ['.crai']}
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  subset_fai: {type: File, doc: "Subsetted fasta index, or any tsv with chrosome names in the first column"}
  snp_vcf: {type: File, doc: "Germline vcf with sites to filter pielup on"}
outputs:
  pileup:
    type: File
    outputBinding:
      glob: $(inputs.input_reads.basename).mini.pileup
