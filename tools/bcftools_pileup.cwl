cwlVersion: v1.0
class: CommandLineTool
id: bcftools_pileup
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

      ${
        var cmd = "ln -s " + inputs.bed_file.path + " pileup.bed";
        if (inputs.read_len != null){
          pad = inputs.read_len/2
          cmd = "cut -f 1,2 " + inputs.reference.secondaryFiles[0].path + " > genome.sizes && bedtools slop -i " + inputs.bed_file.path + " -b " + pad.toString() + " -g genome.sizes | bedtools sort | bedtools merge > pileup.bed"
        }
        return cmd
      }

      bcftools mpileup
      --threads $(inputs.threads)
      -d 8000
      -q 1
      -f $(inputs.reference.path)
      -R pileup.bed
      $(inputs.input_reads.path)
      > $(inputs.input_reads.nameroot).pileup
inputs:
  input_reads: {type: File, secondaryFiles: ['^.bai']}
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  read_len: {type: ['null', int], doc: "set if you want to pad bed file based on read len/2"}
  bed_file: File
outputs:
  pileup:
    type: File
    outputBinding:
      glob: '*.pileup'
