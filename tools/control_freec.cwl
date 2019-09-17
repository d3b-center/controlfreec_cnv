cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-11-5

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'migbro/controlfreec:latest'

baseCommand: [tar, -xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.ref_chr.path)
      && /FREEC/src/freec
      -conf $(inputs.config_file.path)
      && mv $(inputs.tumor_bam.basename)_ratio.txt $(inputs.output_basename).ratio.txt
      && gzip $(inputs.output_basename).ratio.txt
      && mv $(inputs.tumor_bam.basename)_CNVs $(inputs.output_basename).CNVs
      && mv $(inputs.tumor_bam.basename)_sample.cpn $(inputs.output_basename)_sample.cpn
      && gzip $(inputs.output_basename)_sample.cpn
      && mv $(inputs.normal_bam.basename)_control.cpn $(inputs.output_basename)_control.cpn
      && gzip $(inputs.output_basename)_control.cpn

      ${
        var cmd = "mv " + inputs.tumor_bam.basename + "_BAF.txt " + inputs.output_basename + "_tumor_BAF.txt";
        if (inputs.b_allele == null && inputs.tumor_mini_pileup == null){
          cmd = "echo No b allele or mini pileup file, skipping output BAF";
        }
        return cmd;
      }

inputs:
  tumor_bam: { type: File, secondaryFiles: [^.bai] }
  tumor_mini_pileup: {type: ['null', File], doc: "Add if you have a pre-compiled pileup for b allele freq"}
  normal_bam: { type: File, secondaryFiles: [^.bai] }
  normal_mini_pileup: {type: ['null', File], doc: "Add if you have a pre-compiled pileup for b allele freq"}
  ref_chr: {type: File, doc: "folder of reference chromosomes"}
  chr_len: {type: File, doc: "file with chromosome lengths"}
  threads: int
  output_basename: string
  config_file: File
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}
  dbsnp_vcf: {type: ['null', File], , secondaryFiles: [.tbi],  doc: "Gzipped dbsnp for BAF"}
  mappability_file: {type: ['null', File], doc: "GEM output mapability file"}
  reference: {type: ['null', File], doc: "Needed if providing b allele"}
  b_allele: {type: ['null', File], doc: "somatic calls, needed for BAF"}

outputs:
  output_txt:
    type: File
    outputBinding:
      glob: '*.ratio.txt.gz'
  output_cnv:
    type: File
    outputBinding:
      glob: '*.CNVs'
  output_baf:
    type: ['null', File]
    outputBinding:
      glob: '*_tumor_BAF.txt'
  tumor_cpn:
    type: File
    outputBinding:
      glob: '*_sample.cpn.gz'
  normal_cpn:
    type: File
    outputBinding:
      glob: '*_control.cpn.gz'