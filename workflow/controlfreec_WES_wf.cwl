cwlVersion: v1.0
class: Workflow
id: controlfreec_exome_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: { type: File, secondaryFiles: [.crai] }
  input_normal: { type: File, secondaryFiles: [.crai] }
  ref_chr: {type: File, doc: "folder of reference chromosomes"}
  reference: {type: File, secondaryFiles: [.fai]}
  chr_len: {type: File, doc: "file with chromosome lengths"}
  threads: int
  output_basename: string
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}
  exome_flag: {type: string, doc: "insert 'Y' if exome mode"}

outputs:
  output_cnv: {type: File, outputSource: control_free_c/output_cnv}
  output_ratio: {type: File, outputSource: control_free_c/output_txt}

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]
  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]
  gen_config:
    run: ../tools/gen_controlfreec_configfile.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      capture_regions: capture_regions
      exome_flag: exome_flag
      chr_len: chr_len
      threads: threads
    out: [config_file]
  control_free_c: 
    run: ../tools/control_freec.cwl
    in: 
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      capture_regions: capture_regions
      ref_chr: ref_chr
      chr_len: chr_len
      threads: threads
      output_basename: output_basename
      config_file: gen_config/config_file
    out: [output_txt, output_cnv]
    
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2