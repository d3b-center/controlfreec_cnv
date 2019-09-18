cwlVersion: v1.0
class: Workflow
id: controlfreec_bam_config_in_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: {type: File, secondaryFiles: [^.bai]}
  input_normal: {type: File, secondaryFiles: [^.bai]}
  threads: int
  output_basename: string
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_orientation_sample: {type: ['null', string], default: "RF", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  mate_orientation_control: {type: ['null', string], default: "RF", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}
  reference: {type: ['null', File], doc: "Needed if providing b allele"}
  subset_fai: {type: File, doc: "fasta index that is a subset of the main reference fasta file"}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF.  VarDict input recommended"}
  chr_len: {type: File, doc: "TSV with chromsome names and lengths. Limit to chromosome you actualy want analyzed"}

outputs:
  ctrlfreec_cnv: {type: File, outputSource: control_free_c/cnvs}
  ctrlfreec_pval: {type: File, outputSource: control_free_c/cnvs_pvalue}
  ctrlfreec_config: {type: File, outputSource: control_free_c/config_script}
  ctrlfreec_pngs: {type: 'File[]', outputSource: control_free_c/pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: control_free_c/ratio}
  ctrlfreec_baf: {type: File, outputSource: control_free_c/sample_BAF}
  ctrlfreec_info: {type: File, outputSource: control_free_c/info_txt}
  ctrlfreec_tumor_cpn: {type: File, outputSource: control_free_c/sample_cpn}
  ctrlfreec_normal_cpn: {type: File, outputSource: control_free_c/control_cpn}

steps:
  samtools_tumor_pileup:
    run: ../tools/samtools_full_pileup.cwl
    in:
      input_reads: input_tumor
      threads: threads
      reference: reference
      bedtools_genome: subset_fai
    out:
      [pileup]

  samtools_normal_pileup:
    run: ../tools/samtools_full_pileup.cwl
    in:
      input_reads: input_normal
      threads: threads
      reference: reference
      subset_fai: subset_fai
    out:
      [pileup]

  control_free_c: 
    run: ../tools/control-freec-11-6.cwl
    in: 
      mate_file_sample: input_tumor
      mate_orientation_sample: mate_orientation_sample
      sample_pileup: samtools_tumor_pileup/pileup
      mate_file_control: input_normal
      mate_orientation_control: mate_orientation_control
      control_pileup: samtools_normal_pileup/pileup
      chr_len: chr_len
      ploidy: ploidy
      capture_regions: capture_regions
      threads: threads
      reference: reference
      snp_file: b_allele
      output_basename: output_basename
    out: [cnvs, cnvs_pvalue, config_script, pngs, ratio, sample_BAF, info_txt, sample_cpn, control_cpn]
    
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2