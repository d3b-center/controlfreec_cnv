cwlVersion: v1.0
class: Workflow
id: kfdrc_controlfreec_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: {type: File, secondaryFiles: [.crai]}
  input_normal: {type: File, secondaryFiles: [.crai]}
  threads: int
  output_basename: string
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_orientation_sample: {type: ['null', string], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  mate_orientation_control: {type: ['null', string], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}
  reference: {type: File, secondaryFiles: [.fai], doc: "Needed if providing b allele"}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF.  VarDict input recommended.  Tool will prefilter for germline and pass if expression given"}
  chr_len: {type: File, doc: "TSV with chromsome names and lengths. Limit to chromosome you actualy want analyzed"}
  coeff_var: {type: float, default: 0.05, doc: "Coefficient of variantion to set window size.  Default 0.05 recommended"}
  contamination_adjustment: {type: ['null', string], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  include_expression: {type: ['null', string], doc: "Filter expression if vcf has mixed somatic/germline calls, use as-needed"}
  exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has mixed somatic/germline calls, use as-needed"}
  sex: {type: ['null', string], doc: "If known, XX for female, XY for male"}

outputs:
  ctrlfreec_cnv: {type: File, outputSource: rename_outputs/ctrlfreec_cnv}
  ctrlfreec_pval: {type: File, outputSource: rename_outputs/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: rename_outputs/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: rename_outputs/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: rename_outputs/ctrlfreec_bam_ratio}
  ctrlfreec_baf: {type: File, outputSource: rename_outputs/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: rename_outputs/ctrlfreec_info}
  ctrlfreec_tumor_cpn: {type: File, outputSource: rename_outputs/ctrlfreec_tumor_cpn}
  ctrlfreec_normal_cpn: {type: File, outputSource: rename_outputs/ctrlfreec_normal_cpn}

steps:
  bcftools_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: b_allele
      include_expression: include_expression
      exclude_expression: exclude_expression
      output_basename: output_basename
    out:
      [filtered_vcf]

  # samtools_tumor_pileup:
  #   run: ../tools/samtools_full_pileup.cwl
  #   in:
  #     input_reads: input_tumor
  #     threads: threads
  #     reference: reference
  #     subset_fai: subset_fai
  #   out:
  #     [pileup]

  # samtools_normal_pileup:
  #   run: ../tools/samtools_full_pileup.cwl
  #   in:
  #     input_reads: input_normal
  #     threads: threads
  #     reference: reference
  #     subset_fai: subset_fai
  #   out:
  #     [pileup]

  samtools_tumor_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_cram2bam_tumor/bam_file
      threads: threads
      reference: reference
      snp_vcf: b_allele
    out:
      [pileup]

  samtools_normal_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_cram2bam_normal/bam_file
      threads: threads
      reference: reference
      snp_vcf: b_allele
    out:
      [pileup]

  samtools_cram2bam_tumor:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor
      threads: threads
      reference: reference
    out:
      [bam_file]

  samtools_cram2bam_normal:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal
      threads: threads
      reference: reference
    out:
      [bam_file]

  control_free_c: 
    run: ../tools/control-freec-11-6-sbg.cwl
    in: 
      mate_file_sample: samtools_cram2bam_tumor/bam_file
      mate_orientation_sample: mate_orientation_sample
      mini_pileup_sample: samtools_tumor_pileup/pileup
      mate_file_control: samtools_cram2bam_normal/bam_file
      mate_orientation_control: mate_orientation_control
      mini_pileup_control: samtools_normal_pileup/pileup
      chr_len: chr_len
      ploidy: ploidy
      capture_regions: capture_regions
      max_threads: threads
      reference: reference
      snp_file: bcftools_filter_vcf/filtered_vcf
      coeff_var: coeff_var
      output_basename: output_basename
      sex: sex
      contamination_adjustment: contamination_adjustment
    out: [cnvs, cnvs_pvalue, config_script, pngs, ratio, sample_BAF, info_txt, sample_cpn, control_cpn]

  rename_outputs:
    run: ../tools/ubuntu_rename_outputs.cwl
    in:
      input_files: [control_free_c/cnvs, control_free_c/cnvs_pvalue, control_free_c/config_script, control_free_c/ratio, control_free_c/sample_BAF, control_free_c/info_txt, control_free_c/sample_cpn, control_free_c/control_cpn]
      input_pngs: control_free_c/pngs
    out: [ctrlfreec_cnv, ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_baf, ctrlfreec_info, ctrlfreec_tumor_cpn, ctrlfreec_normal_cpn]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4