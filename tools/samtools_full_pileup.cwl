cwlVersion: v1.0
class: CommandLineTool
id: samtools_full_pileup
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 16
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      cut -f 1  $(inputs.subset_fai.path)
      | xargs -IXX -P 16 sh -c  "samtools mpileup -f $(inputs.reference.path) -r XX -d 250 $(inputs.input_reads.path)
      | awk '{if(\$4 != 0) print \$0}' > XX.pileup"

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
    value: c5.4xlarge;ebs-gp2;2000
