cwlVersion: v1.0
class: Workflow
id: controlfreec_bam_config_in_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: { type: File, secondaryFiles: [^.bai] }
  input_normal: { type: File, secondaryFiles: [^.bai] }
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
  output_cnv: {type: File, outputSource: control_free_c/output_cnv}
  ctrlfreec_bam_ratio: {type: File, outputSource: control_free_c/output_txt}
  ctrlfreec_pval: {type: File, outputSource: control_free_c_r/output_pval}
  ctrlfreec_png: {type: File, outputSource: control_free_c_viz/output_png}
  ctrlfreec_baf: {type: File, outputSource: control_free_c/output_baf}
  ctrlfreec_tumor_cpn: {type: File, outputSource: control_free_c/tumor_cpn}
  ctrlfreec_normal_cpn: {type: File, outputSource: control_free_c/normal_cpn}

steps:
  control_free_c: 
    run: ../tools/control_freec.cwl
    in: 
      tumor_bam: input_tumor
      normal_bam: input_normal
      capture_regions: capture_regions
      ref_chr: ref_chr
      chr_len: chr_len
      threads: threads
      output_basename: output_basename
      config_file: config_file
      dbsnp_vcf: dbsnp_vcf
      mappability_file: mappability_file
      reference: reference
      b_allele: b_allele
    out: [output_txt, output_cnv, output_baf, tumor_cpn, normal_cpn]
  control_free_c_r:
    run: ../tools/control_freec_R.cwl
    in:
      cnv_bam_ratio: control_free_c/output_txt
      cnv_result: control_free_c/output_cnv
    out: [output_pval]
  control_free_c_viz:
    run: ../tools/control_freec_visualize.cwl
    in:
      output_basename: output_basename
      cnv_bam_ratio: control_free_c/output_txt
    out: [output_png]
    
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2