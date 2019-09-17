cwlVersion: v1.0
class: CommandLineTool
id: samtools_full_pileup
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

      cut -f 1  $(inputs.bedtools_genome.path) 
      | xargs -IXX -P $(inputs.threads) sh -c  "samtools mpileup -f $(inputs.reference.path) -r XX -d 8000 -q 1 -Q 0 $(inputs.input_reads.path)
      > XX.pileup"

      cut -f 1 $(inputs.bedtools_genome.path) | xargs -IXX cat XX.pileup >> $(inputs.input_reads.basename).pileup

inputs:
  input_reads: {type: File, secondaryFiles: ['^.bai']}
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  bedtools_genome: {type: File, doc: "File with chromosomes and sizes of bed file contents"}
outputs:
  pileup:
    type: File
    outputBinding:
      glob: $(inputs.input_reads.basename).pileup
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge;ebs-gp2;2000
