cwlVersion: v1.0
class: CommandLineTool
id: samtools_pileup
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

      bedtools sort -i $(inputs.bed_file.path) -faidx $(inputs.bedtools_genome.path) > $(inputs.bed_file.nameroot).safety_sorted.bed

      cut -f 1  $(inputs.bedtools_genome.path) 
      | xargs -IXX -P $(inputs.threads) sh -c  "samtools mpileup -d 8000 -q 1 -f $(inputs.reference.path) -r XX -Q 0 $(inputs.input_reads.path)
      | awk {'printf (\"%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", \$1,\$2-1,\$2,\$3,\$4,\$5,\$6)'}
      | bedtools intersect -a stdin -b $(inputs.bed_file.nameroot).safety_sorted.bed -g $(inputs.bedtools_genome.path) -sorted -wa
      | cut -f 1,3- > XX.pileup"

      cut -f 1 $(inputs.bedtools_genome.path) | xargs -IXX cat XX.pileup >> $(inputs.input_reads.basename).mini.pileup

inputs:
  input_reads: {type: File, secondaryFiles: ['^.bai']}
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  bedtools_genome: {type: File, doc: "File with chromosomes and sizes of bed file contents"}
  bed_file: File
outputs:
  pileup:
    type: File
    outputBinding:
      glob: $(inputs.input_reads.basename).mini.pileup
