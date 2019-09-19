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
    coresMin: 36
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      cut -f 1  $(inputs.subset_fai.path) 
      | xargs -IXX -P 36 sh -c  "samtools mpileup -f $(inputs.reference.path) -r XX -d 8000 -q 1 -Q 0 $(inputs.input_reads.path)
      > XX.pileup"

      cut -f 1 $(inputs.subset_fai.path) | xargs -IXX cat XX.pileup >> $(inputs.input_reads.basename).pileup

inputs:
  input_reads: {type: File, secondaryFiles: ['.crai']}

  reference: {type: File, secondaryFiles: [.fai]}
  subset_fai: {type: File, doc: "Subset of reference fai with desired chromosomes to split out"}
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
