class: Workflow
cwlVersion: v1.0
id: brownm28/mb-controlfreec-troubleshoot/control-freec-loh-11-6/4
label: Control-FREEC LOH 11.6
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: filtered_reference_fai
    'sbg:fileTypes': FAI
    type: File?
    label: Filtered Reference FAI
    'sbg:x': -599.25
    'sbg:y': -225
  - id: tumor_bam
    'sbg:fileTypes': BAM
    type: File
    label: Tumor BAM
    'sbg:x': -714.75
    'sbg:y': -366
  - id: normal_bam
    'sbg:fileTypes': BAM
    type: File
    label: Normal BAM
    'sbg:x': -712.75
    'sbg:y': -84
  - id: reference_fasta
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: Reference FASTA
    'sbg:x': -596
    'sbg:y': -475
  - id: chromosome_sizes
    'sbg:fileTypes': 'TXT, LEN, SIZES'
    type: File
    label: Chromosome Sizes
    'sbg:x': -570
    'sbg:y': 19
  - id: capture_regions
    'sbg:fileTypes': BED
    type: File?
    label: Capture Regions
    'sbg:x': -430
    'sbg:y': 90
  - id: snp_file
    'sbg:fileTypes': 'TXT, VCF'
    type: File?
    label: SNP File
    'sbg:x': -464
    'sbg:y': -545
  - id: mate_orientation_control
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - RF
          - FR
          - FF
        name: mate_orientation_control
    'sbg:exposed': true
  - id: mate_orientation_sample
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - RF
          - FR
          - FF
        name: mate_orientation_sample
    'sbg:exposed': true
  - id: noisy_data
    type: boolean?
    'sbg:exposed': true
  - id: ploidy
    type: 'int[]'
    'sbg:exposed': true
  - id: sex
    type:
      - 'null'
      - type: enum
        symbols:
          - XX
          - XY
        name: sex
    'sbg:exposed': true
outputs:
  - id: cnvs_pvalue
    outputSource:
      - control_freec_11_5/cnvs_pvalue
    'sbg:fileTypes': TXT
    type: File?
    'sbg:x': 259.1755676269531
    'sbg:y': -284.083984375
  - id: config_script
    outputSource:
      - control_freec_11_5/config_script
    'sbg:fileTypes': TXT
    type: File?
    'sbg:x': 226
    'sbg:y': -394
  - id: sample_BAF
    outputSource:
      - control_freec_11_5/sample_BAF
    'sbg:fileTypes': TXT
    type: File?
    'sbg:x': 195.3478240966797
    'sbg:y': -514.7825927734375
  - id: ratio
    outputSource:
      - control_freec_11_5/ratio
    'sbg:fileTypes': TXT
    type: File?
    'sbg:x': 257.6692810058594
    'sbg:y': -122.2210464477539
  - id: sample_pileup
    outputSource:
      - control_freec_11_5/sample_pileup
    'sbg:fileTypes': PILEUP
    type: File?
    'sbg:x': -37.18376922607422
    'sbg:y': -507.1328430175781
  - id: sample_cpn
    outputSource:
      - control_freec_11_5/sample_cpn
    'sbg:fileTypes': CPN
    type: File?
    'sbg:x': 86.4817123413086
    'sbg:y': -685.0687866210938
  - id: pngs
    outputSource:
      - control_freec_11_5/pngs
    'sbg:fileTypes': PNG
    type: 'File[]?'
    'sbg:x': 459.25750732421875
    'sbg:y': -232.2218017578125
  - id: control_pileup
    outputSource:
      - control_freec_11_5/control_pileup
    'sbg:fileTypes': PILEUP
    type: File?
    'sbg:x': 534.8803100585938
    'sbg:y': -102.32856750488281
  - id: cnvs
    outputSource:
      - control_freec_11_5/cnvs
    'sbg:fileTypes': TXT
    type: File?
    'sbg:x': 186.12583923339844
    'sbg:y': 35.57179260253906
  - id: control_cpn
    outputSource:
      - control_freec_11_5/control_cpn
    'sbg:fileTypes': CPN
    type: File?
    'sbg:x': 498.4034118652344
    'sbg:y': 92.51129150390625
steps:
  - id: samtools_mpileup_parallel_1_6_tumor
    in:
      - id: fai_file
        source: filtered_reference_fai
      - id: input_bam_files
        source:
          - tumor_bam
      - id: reference_fasta
        source: reference_fasta
    out:
      - id: output_file
    run:
      class: Workflow
      cwlVersion: v1.0
      id: vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3
      label: SAMtools Mpileup Parallel 1.6
      $namespaces:
        sbg: 'https://sevenbridges.com'
      inputs:
        - id: fai_file
          'sbg:fileTypes': FAI
          type: File?
          label: Reference FAI
          'sbg:x': -509
          'sbg:y': -119
        - id: input_bam_files
          'sbg:fileTypes': BAM
          type: 'File[]'
          label: BAM files
          'sbg:x': -512
          'sbg:y': 150
        - id: reference_fasta
          'sbg:fileTypes': 'FASTA, FA, GZ'
          type: File?
          label: Reference FASTA
          'sbg:x': -508
          'sbg:y': 18
      outputs:
        - id: output_file
          outputSource:
            - sbg_samtools_merge_mpileup/output_file
          'sbg:fileTypes': 'PILEUP, BCF, VCF'
          type: File?
          label: Output file
          'sbg:x': 373.5552978515625
          'sbg:y': 27.004995346069336
      steps:
        - id: sbg_prepare_intervals
          in:
            - id: fai_file
              source: fai_file
            - id: format
              default: chr start end
            - id: split_mode
              default: File per chr with alt contig in a single file
          out:
            - id: intervals
            - id: names
            - id: str_arr
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/sbg-prepare-intervals/0
            baseCommand:
              - python
              - sbg_prepare_intervals.py
            inputs:
              - 'sbg:category': File Input
                id: bed_file
                type: File?
                inputBinding:
                  position: 1
                  prefix: '--bed'
                  shellQuote: false
                label: Input BED file
                doc: >-
                  Input BED file containing intervals. Required for modes 3 and
                  4.
                'sbg:fileTypes': BED
              - 'sbg:category': File Input
                id: fai_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '--fai'
                  shellQuote: false
                label: Input FAI file
                doc: >-
                  FAI file is converted to BED format if BED file is not
                  provided.
                'sbg:fileTypes': FAI
              - 'sbg:category': Input
                id: format
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - chr start end
                      - 'chr:start-end'
                    name: format
                label: Interval format
                doc: Format of the intervals in the generated files.
              - default: 0
                'sbg:category': Input
                id: split_mode
                type:
                  type: enum
                  symbols:
                    - File per interval
                    - File per chr with alt contig in a single file
                    - Output original BED
                    - File per interval with alt contig in a single file
                  name: split_mode
                inputBinding:
                  position: 3
                  prefix: '--mode'
                  shellQuote: false
                  valueFrom: |-
                    ${
                        if (self == 0) {
                            self = null;
                            inputs.split_mode = null
                        };


                        mode = inputs.split_mode
                        switch (mode) {
                            case "File per interval":
                                return 1
                            case "File per chr with alt contig in a single file":
                                return 2
                            case "Output original BED":
                                return 3
                            case "File per interval with alt contig in a single file":
                                return 4
                        }
                        return 3
                    }
                label: Split mode
                doc: >-
                  Depending on selected Split Mode value, output files are
                  generated in accordance with description below:  1. File per
                  interval - The tool creates one interval file per line of the
                  input BED(FAI) file. Each interval file contains a single line
                  (one of the lines of BED(FAI) input file).  2. File per chr
                  with alt contig in a single file - For each contig(chromosome)
                  a single file is created containing all the intervals
                  corresponding to it . All the intervals (lines) other than
                  (chr1, chr2 ... chrY or 1, 2 ... Y) are saved as
                  ("others.bed").  3. Output original BED - BED file is required
                  for execution of this mode. If mode 3 is applied input is
                  passed to the output.  4. File per interval with alt contig in
                  a single file - For each chromosome a single file is created
                  for each interval. All the intervals (lines) other than (chr1,
                  chr2 ... chrY or 1, 2 ... Y) are saved as ("others.bed").
                  NOTE: Do not use option 1 (File per interval) with exome BED
                  or a BED with a lot of GL contigs, as it will create a large
                  number of files.
            outputs:
              - id: intervals
                doc: Array of BED files genereted as per selected Split Mode.
                label: Intervals
                type: 'File[]?'
                outputBinding:
                  glob: Intervals/*.bed
                  outputEval: |-
                    ${

                        for (var i = 0; i < self.length; i++) {
                            var out_metadata = {
                                'sbg_scatter': 'true'
                            };
                            self[i] = setMetadata(self[i], out_metadata)
                        };

                        return self

                    }
                'sbg:fileTypes': BED
              - id: names
                doc: File containing the names of created files.
                label: Output file names
                type: string?
                outputBinding:
                  loadContents: true
                  glob: Intervals/names.txt
                  outputEval: |-
                    ${
                        content = self[0].contents.replace(/\0/g, '')
                        content = content.replace('[', '')
                        content = content.replace(']', '')
                        content = content.replace(/\'/g, "")
                        content = content.replace(/\s/g, '')
                        content_arr = content.split(",")

                        return content_arr


                    }
              - id: str_arr
                doc: Outputs BED content as strings.
                label: String output
                type: 'string[]?'
                outputBinding:
                  loadContents: true
                  glob: |-
                    ${
                        if (inputs.bed_file) {
                            glob = inputs.bed_file.path
                            glob = glob.split('/').slice(-1)[0]
                        } else if (inputs.fai_file) {
                            glob = inputs.fai_file.path
                            glob = glob.split('/').slice(-1)[0].split('.').slice(0, -1).join('.') + '.bed'
                        }

                        return glob
                    }
                  outputEval: |-
                    ${
                        rows = self[0].contents
                        if (rows[rows.length - 1] == '\n') {
                            rows = rows.split(/\r?\n/).slice(0, -1);
                        } else {
                            rows = rows.split(/\r?\n/);
                        }
                        out_list = []
                        for (i = 0; i < rows.length; i++) {
                            row = rows[i];
                            chromosome = row.split("\t")[0];
                            start = row.split("\t")[1];
                            end = row.split("\t")[2];
                            if (start) {
                                interval = chromosome.concat(":", start, "-", end);
                            } else {
                                interval = chromosome
                            }
                            out_list.push(interval);
                        }
                        return out_list;

                    }
            doc: >-
              Depending on selected Split Mode value, output files are generated
              in accordance with description below:


              1. File per interval - The tool creates one interval file per line
              of the input BED(FAI) file.

              Each interval file contains a single line (one of the lines of
              BED(FAI) input file).


              2. File per chr with alt contig in a single file - For each
              contig(chromosome) a single file

              is created containing all the intervals corresponding to it .

              All the intervals (lines) other than (chr1, chr2 ... chrY or 1, 2
              ... Y) are saved as

              ("others.bed").


              3. Output original BED - BED file is required for execution of
              this mode. If mode 3 is applied input is passed to the output.


              4. File per interval with alt contig in a single file - For each
              chromosome a single file is created for each interval.

              All the intervals (lines) other than (chr1, chr2 ... chrY or 1, 2
              ... Y) are saved as

              ("others.bed").


              ##### Common issues: 

              Do not use option 1 (File per interval) with exome BED or a BED
              with a lot of GL contigs, as it will create a large number of
              files.
            label: SBG Prepare Intervals
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: |-
                  ${
                      if (inputs.format)
                          return "--format " + "\"" + inputs.format + "\""
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'images.sbgenomics.com/bogdang/sbg_prepare_intervals:1.0'
              - class: InitialWorkDirRequirement
                listing:
                  - entryname: sbg_prepare_intervals.py
                    entry: >-
                      """

                      Usage:
                          sbg_prepare_intervals.py [options] [--fastq FILE --bed FILE --mode INT --format STR --others STR]

                      Description:
                          Purpose of this tool is to split BED file into files based on the selected mode.
                          If bed file is not provided fai(fasta index) file is converted to bed.

                      Options:

                          -h, --help            Show this message.

                          -v, -V, --version     Tool version.

                          -b, -B, --bed FILE    Path to input bed file.

                          --fai FILE            Path to input fai file.

                          --format STR          Output file format.

                          --mode INT            Select input mode.

                      """



                      import os

                      import sys

                      import glob

                      import shutil

                      from docopt import docopt


                      default_extension = '.bed'  # for output files



                      def create_file(contents, contig_name,
                      extension=default_extension):
                          """function for creating a file for all intervals in a contig"""

                          new_file = open("Intervals/" + contig_name + extension, "w")
                          new_file.write(contents)
                          new_file.close()


                      def add_to_file(line, name, extension=default_extension):
                          """function for adding a line to a file"""

                          new_file = open("Intervals/" + name + extension, "a")
                          if lformat == formats[1]:
                              sep = line.split("\t")
                              line = sep[0] + ":" + sep[1] + "-" + sep[2]
                          new_file.write(line)
                          new_file.close()


                      def fai2bed(fai):
                          """function to create a bed file from fai file"""

                          region_thr = 10000000  # threshold used to determine starting point accounting for telomeres in chromosomes
                          basename = fai[0:fai.rfind(".")]
                          with open(fai, "r") as ins:
                              new_array = []
                              for line in ins:
                                  len_reg = int(line.split()[1])
                                  cutoff = 0 if (
                                  len_reg < region_thr) else 0  # sd\\telomeres or start with 1
                                  new_line = line.split()[0] + '\t' + str(cutoff) + '\t' + str(
                                      len_reg + cutoff)
                                  new_array.append(new_line)
                          new_file = open(basename + ".bed", "w")
                          new_file.write("\n".join(new_array))
                          return basename + ".bed"


                      def chr_intervals(no_of_chrms=23):
                          """returns all possible designations for chromosome intervals"""

                          chrms = []
                          for i in range(1, no_of_chrms):
                              chrms.append("chr" + str(i))
                              chrms.append(str(i))
                          chrms.extend(["x", "y", "chrx", "chry"])
                          return chrms


                      def mode_1(orig_file):
                          """mode 1: every line is a new file"""

                          with open(orig_file, "r") as ins:
                              prev = ""
                              counter = 0
                              names = []
                              for line in ins:
                                  if is_header(line):
                                      continue
                                  if line.split()[0] == prev:
                                      counter += 1
                                  else:
                                      counter = 0
                                  suffix = "" if (counter == 0) else "_" + str(counter)
                                  create_file(line, line.split()[0] + suffix)
                                  names.append(line.split()[0] + suffix)
                                  prev = line.split()[0]

                              create_file(str(names), "names", extension=".txt")


                      def mode_2(orig_file, others_name):
                          """mode 2: separate file is created for each chromosome, and one file is created for other intervals"""

                          chrms = chr_intervals()
                          names = []

                          with open(orig_file, 'r') as ins:
                              for line in ins:
                                  if is_header(line):
                                      continue
                                  name = line.split()[0]
                                  if name.lower() in chrms:
                                      name = name
                                  else:
                                      name = others_name
                                  try:
                                      add_to_file(line, name)
                                      if not name in names:
                                          names.append(name)
                                  except:
                                      raise Exception(
                                          "Couldn't create or write in the file in mode 2")

                              create_file(str(names), "names", extension=".txt")


                      def mode_3(orig_file, extension=default_extension):
                          """mode 3: input file is staged to output"""

                          orig_name = orig_file.split("/")[len(orig_file.split("/")) - 1]
                          output_file = r"./Intervals/" + orig_name[
                                                          0:orig_name.rfind('.')] + extension

                          shutil.copyfile(orig_file, output_file)

                          names = [orig_name[0:orig_name.rfind('.')]]
                          create_file(str(names), "names", extension=".txt")


                      def mode_4(orig_file, others_name):
                          """mode 4: every interval in chromosomes is in a separate file. Other intervals are in a single file"""

                          chrms = chr_intervals()
                          names = []

                          with open(orig_file, "r") as ins:
                              counter = {}
                              for line in ins:
                                  if line.startswith('@'):
                                      continue
                                  name = line.split()[0].lower()
                                  if name in chrms:
                                      if name in counter:
                                          counter[name] += 1
                                      else:
                                          counter[name] = 0
                                      suffix = "" if (counter[name] == 0) else "_" + str(counter[name])
                                      create_file(line, name + suffix)
                                      names.append(name + suffix)
                                      prev = name
                                  else:
                                      name = others_name
                                      if not name in names:
                                          names.append(name)
                                      try:
                                          add_to_file(line, name)
                                      except:
                                          raise Exception(
                                              "Couldn't create or write in the file in mode 4")

                          create_file(str(names), "names", extension=".txt")


                      def prepare_intervals():
                          # reading input files and split mode from command line
                          args = docopt(__doc__, version='1.0')

                          bed_file = args['--bed']
                          fai_file = args['--fai']
                          split_mode = int(args['--mode'])

                          # define file name for non-chromosomal contigs
                          others_name = 'others'

                          global formats, lformat
                          formats = ["chr start end", "chr:start-end"]
                          lformat = args['--format']
                          if lformat == None:
                              lformat = formats[0]
                          if not lformat in formats:
                              raise Exception('Unsuported interval format')

                          if not os.path.exists(r"./Intervals"):
                              os.mkdir(r"./Intervals")
                          else:
                              files = glob.glob(r"./Intervals/*")
                              for f in files:
                                  os.remove(f)

                          # create variable input_file taking bed_file as priority
                          if bed_file:
                              input_file = bed_file
                          elif fai_file:
                              input_file = fai2bed(fai_file)
                          else:
                              raise Exception('No input files are provided')

                          # calling adequate split mode function
                          if split_mode == 1:
                              mode_1(input_file)
                          elif split_mode == 2:
                              mode_2(input_file, others_name)
                          elif split_mode == 3:
                              if bed_file:
                                  mode_3(input_file)
                              else:
                                  raise Exception('Bed file is required for mode 3')
                          elif split_mode == 4:
                              mode_4(input_file, others_name)
                          else:
                              raise Exception('Split mode value is not set')


                      def is_header(line):
                          x = line.split('\t')
                          try:
                              int(x[1])
                              int(x[2])
                              header = False
                          except:
                              sys.stderr.write('Line is skipped: {}'.format(line))
                              header = True
                          return header


                      if __name__ == '__main__':
                          prepare_intervals()
                    writable: false
                  - $(inputs.fai_file)
                  - $(inputs.bed_file)
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - Converters
            'sbg:cmdPreview': python sbg_prepare_intervals.py  --format "chr start end" --mode 2
            'sbg:image_url': null
            'sbg:license': Apache License 2.0
            'sbg:toolAuthor': Seven Bridges Genomics
            'sbg:toolkit': SBGTools
            'sbg:toolkitVersion': '1.0'
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/sbg-prepare-intervals/0
            'sbg:revision': 0
            'sbg:revisionNotes': >-
              Copy of
              vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
            'sbg:modifiedOn': 1560510063
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510063
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 0
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510063
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
            'sbg:publisher': sbg
            'sbg:content_hash': aa51632b14cdb5f6e798ed50dbb067444bf78dfa99bbdb384445c93e9f8782262
            'sbg:copyOf': >-
              vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
          label: SBG Prepare Intervals
          'sbg:x': -294
          'sbg:y': -63
        - id: samtools_mpileup_1_6
          in:
            - id: input_bam_files
              source:
                - input_bam_files
            - id: output_format
              default: PILEUP
            - id: reference_fasta_file
              source: reference_fasta
            - id: region_pileup_generation
              source: sbg_prepare_intervals/names
          out:
            - id: output_pileup_vcf_or_bcf_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/samtools-mpileup-1-6/1
            baseCommand: []
            inputs:
              - 'sbg:category': Input options
                id: assume_the_quality_Illumina_encoding
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-6'
                  shellQuote: false
                label: Assume the quality is in the Illumina-1.3+ encoding
                doc: Assume the quality is in the Illumina-1.3+ encoding.
              - 'sbg:category': File input
                id: chr_pos_list_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-l'
                  shellQuote: false
                label: List of positions (chr pos) or regions in BED file
                doc: >-
                  BED or position list file containing a list of regions or
                  sites where pileup or BCF should be generated.
                'sbg:fileTypes': BED
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: comma_separated_list_of_platforms_for_indels
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-P'
                  shellQuote: false
                label: Comma separated list of platforms for indels
                doc: >-
                  Comma dilimited list of platforms (determined by @RG-PL) from
                  which indel candidates are obtained. It is recommended to
                  collect indel candidates from sequencing technologies that
                  have low indel error rate such as ILLUMINA.
              - 'sbg:category': Input options
                id: count_anomalous_read_pairs
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-A'
                  shellQuote: false
                label: Count anomalous read pairs
                doc: Do not skip anomalous read pairs in variant calling.
              - 'sbg:category': Input options
                id: disable_baq_computation
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-B'
                  shellQuote: false
                label: Disable BAQ computation
                doc: >-
                  Disable probabilistic realignment for the computation of base
                  alignment quality (BAQ). BAQ is the Phred-scaled probability
                  of a read base being misaligned. Applying this option greatly
                  helps to reduce false SNPs caused by misalignments.
              - 'sbg:category': File input
                id: exclude_read_groups_list_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-G'
                  shellQuote: false
                label: Exclude read groups listed in file
                doc: Exclude read groups listed in file.
                'sbg:fileTypes': TXT
              - 'sbg:category': Input options
                id: filter_flags
                type: string?
                inputBinding:
                  position: 2
                  prefix: '--ff'
                  shellQuote: false
                label: Filter flags
                doc: >-
                  Filter flags: skip reads with mask bits set
                  [UNMAP,SECONDARY,QCFAIL,DUP].
              - 'sbg:category': output configuration
                id: filter_zero_coverage_lines
                type: boolean?
                label: Filter zero coverage lines.
                doc: >-
                  Filter zero coverage lines from pileup files. This option is
                  valid only when output pileup files.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: gap_extension_seq_error_probability
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-e'
                  shellQuote: false
                label: Phred-scaled gap extension seq error probability
                doc: >-
                  Phred-scaled gap extension sequencing error probability.
                  Reducing INT leads to longer indels.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: gap_open_sequencing_error_probability
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-o'
                  shellQuote: false
                label: Phred-scaled gap open sequencing error probability
                doc: >-
                  Phred-scaled gap open sequencing error probability. Reducing
                  INT leads to more indel calls.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: homopolymer_err_coeficient
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-h'
                  shellQuote: false
                label: Coefficient for homopolymer errors
                doc: >-
                  Coefficient for modeling homopolymer errors. Given an l-long
                  homopolymer run, the sequencing error of an indel of size s is
                  modeled as INT*s/l.
              - 'sbg:category': configuration
                id: ignore_overlaps
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-x'
                  shellQuote: false
                label: Disable read-pair overlap detection
                doc: Disable read-pair overlap detection.
              - 'sbg:category': Input options
                id: ignore_rg_tags
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-R'
                  shellQuote: false
                label: Ignore RG tags
                doc: Ignore rg tags.
              - 'sbg:category': Input options
                schema:
                  - File
                id: input_bam_files
                type: 'File[]'
                inputBinding:
                  position: 3
                  shellQuote: false
                  separator: ' '
                  streamable: false
                label: List of input BAM files
                doc: List of input BAM files.
                'sbg:fileTypes': BAM
                secondaryFiles:
                  - ^.bai
              - 'sbg:category': Input options
                id: mapq_threshold
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-q'
                  shellQuote: false
                  separator: ' '
                label: Skip alignments with mapQ  smaller than
                doc: Minimum mapping quality for an alignment to be used.
              - 'sbg:category': configuration
                schema:
                  - 'null'
                  - int
                id: max_idepth
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-L'
                  shellQuote: false
                label: Skip INDEL calling if the average per-sample depth is above
                doc: Skip INDEL calling if the average per-sample depth is above.
              - 'sbg:category': Input options
                id: max_per_bam_depth
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-d'
                  shellQuote: false
                label: Max per-BAM depth
                doc: 'At a position, read maximally INT reads per input BAM.'
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: min_gapped_reads_for_indel
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-m'
                  shellQuote: false
                label: Minimum gapped reads for indel candidates
                doc: Minimum gapped reads for indel candidates.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: minimum_fraction_of_gapped_reads
                type: float?
                inputBinding:
                  position: 2
                  prefix: '-F'
                  shellQuote: false
                label: Minimum fraction of gapped reads for candidates
                doc: Minimum fraction of gapped reads for candidates.
              - 'sbg:category': configuration
                id: more_info_to_output
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-t'
                  shellQuote: false
                label: Comma-separated list of FORMAT and INFO tags to output
                doc: >-
                  Comma-separated list of FORMAT and INFO tags
                  (DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []") to output.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: no_indel_calling
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-I'
                  shellQuote: false
                label: Do not perform indel calling
                doc: Do not perform indel calling.
              - 'sbg:category': Output options
                id: output_base_positions_on_reads
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-O'
                  shellQuote: false
                label: Output base positions on reads
                doc: Output base positions on reads (disabled by -g/-u).
              - 'sbg:category': configuration
                id: output_format
                type:
                  type: enum
                  symbols:
                    - PILEUP
                    - BCF
                    - VCF
                  name: output_format
                label: Output file format
                doc: Output file format.
              - 'sbg:category': Output options
                id: output_mapping_quality
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-s'
                  shellQuote: false
                label: Output mapping quality
                doc: Output mapping quality (disabled by -g/-u).
              - 'sbg:category': Output options
                id: output_uncompressed_bcf_or_vcf
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-u'
                  shellQuote: false
                label: Generate uncompress BCF/VCF output
                doc: >-
                  Similar to -g except that the output is uncompressed BCF,
                  which is preferred for piping.
              - 'sbg:category': Input options
                id: parameter_for_adjusting_mapq
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-C'
                  shellQuote: false
                label: Parameter for adjusting mapQ
                doc: >-
                  Coefficient for downgrading mapping quality for reads
                  containing excessive mismatches. Given a read with a
                  phred-scaled probability q of being generated from the mapped
                  position, the new mapping quality is about
                  sqrt((INT-q)/INT)*INT. A zero value disables this
                  functionality; if enabled, the recommended value for BWA is
                  50.
              - 'sbg:category': configuration
                id: per_sample_mF
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-p'
                  shellQuote: false
                label: Apply -m and -F thresholds per sample
                doc: >-
                  Apply -m and -F thresholds per sample to increase sensitivity
                  of calling. By default both options are applied to reads
                  pooled from all samples.
              - 'sbg:category': Input options
                id: recalc_BAQ_fly
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-E'
                  shellQuote: false
                  separator: ' '
                label: Recalculate BAQ on the fly
                doc: 'Recalculate BAQ on the fly, ignore existing BQ tags.'
              - 'sbg:category': File input
                schema:
                  - 'null'
                  - File
                id: reference_fasta_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-f'
                  shellQuote: false
                label: Reference FASTA file
                doc: Faidx indexed reference sequence (FASTA) file.
                'sbg:fileTypes': 'FASTA, FA, GZ'
                secondaryFiles:
                  - .fai
              - 'sbg:category': Input options
                'sbg:includeInPorts': true
                id: region_pileup_generation
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-r'
                  shellQuote: false
                label: Region in which pileup is generated
                doc: Region in which pileup is generated.
              - 'sbg:category': Input options
                id: required_flags
                type: int?
                inputBinding:
                  position: 2
                  prefix: '--rf'
                  shellQuote: false
                label: Required flags
                doc: 'Required flags: skip reads with mask bits unset.'
              - 'sbg:category': Input options
                id: skip_bases_with_baq_smaller_than_defined
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-Q'
                  shellQuote: false
                  separator: ' '
                label: Skip bases with baseQ/BAQ smaller than
                doc: Minimum base quality for a base to be considered.
            outputs:
              - id: output_pileup_vcf_or_bcf_file
                doc: 'Output PILEUP, VCF, or BCF file.'
                label: 'Output PILEUP, VCF, or BCF file'
                type: File?
                outputBinding:
                  glob: '{*.pileup,*.bcf,*.vcf}'
                  outputEval: |-
                    ${
                        self = inheritMetadata(self, inputs.input_bam_files);


                        var add_metadata_key_ScatteredUsing = function(self, inputs) {
                            if (inputs.chr_pos_list_file) {
                                len = inputs.chr_pos_list_file.path.split('/').length

                                return inputs.chr_pos_list_file.path.split('/')[len - 1].split('.')[0]
                            } else if (inputs.region_pileup_generation) {
                                return inputs.region_pileup_generation
                            } else {
                                return
                            }
                        }
                        for (var i = 0; i < self.length; i++) {
                            var out_metadata = {
                                'ScatteredUsing': add_metadata_key_ScatteredUsing(self[i], inputs)
                            };
                            self[i] = setMetadata(self[i], out_metadata)
                        };

                        return self

                    }
                'sbg:fileTypes': 'PILEUP, BCF, VCF'
            doc: >-
              SAMtools Mpileup generates BCF or PILEUP for one or multiple BAM
              files. Alignment records are grouped by sample identifiers in @RG
              header lines. If sample identifiers are absent, each input file is
              regarded as one sample.


              In the pileup format (without -uor-g), each line represents a
              genomic position consisting of chromosome name, coordinate,
              reference base, read bases, read qualities, and alignment mapping
              qualities. Information on match, mismatch, indel, strand, mapping
              quality, and the start and end of a read are all encoded at the
              read base column. In this column, a dot stands for a match to the
              reference base on the forward strand, a comma for a match on the
              reverse strand, a '>' or '<' for a reference skip, 'ACGTN' for a
              mismatch on the forward strand, and 'acgtn' for a mismatch on the
              reverse strand. A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there
              is an insertion between this reference position and the next
              reference position. The length of the insertion is given by the
              integer in the pattern followed by the inserted sequence.
              Similarly, a pattern '-[0-9]+[ACGTNacgtn]+' represents a deletion
              from the reference. The deleted bases will be presented as `*' in
              the following lines. Also at the read base column, a symbol '^',
              marks the start of a read. The ASCII of the character following
              '^' minus 33 gives the mapping quality. The symbol '$' marks the
              end of a read segment.


              #### Common Issue:

              - Please use the public pileup parallel pipeline for large input
              files.

              - Please sort BAM files by coordinates before using mpileup.
            label: SAMtools Mpileup 1.6
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: /opt/samtools-1.6/samtools
              - position: 1
                shellQuote: false
                valueFrom: mpileup
              - position: 2
                shellQuote: false
                valueFrom: |-
                  ${

                      if (inputs.output_format) {
                          if (inputs.output_format == "PILEUP")
                              return
                          else if (inputs.output_format == "VCF")
                              return "-v"
                          else if (inputs.output_format == "BCF")
                              return "-g"
                          else
                              return


                      }



                  }
              - position: 1002
                shellQuote: false
                valueFrom: |-
                  ${

                      // Check if input files are delivered in an array or not
                      if (inputs.input_bam_files.constructor == Array) {
                          // Select the first entry of the array
                          filepath = "/" + inputs.input_bam_files[0].path
                      } else

                      {
                          // Use the file name directly
                          filepath = "/" + inputs.input_bam_files.path
                      }


                      filename = filepath.split("/").pop()
                      new_filename = filename.substr(0, filename.lastIndexOf("."))


                      if (inputs.chr_pos_list_file)
                          new_filename = new_filename + '_' + inputs.chr_pos_list_file.path.split('/')[inputs.chr_pos_list_file.path.split('/').length - 1].split('.')[0]

                      if (inputs.region_pileup_generation)
                          new_filename = new_filename + '_' + inputs.region_pileup_generation

                      extension = '.pileup'

                      if ((inputs.output_format && inputs.output_format === "BCF"))
                          extension = '.bcf'

                      else if ((inputs.output_format && inputs.output_format === "VCF"))
                          extension = '.vcf'

                      if (inputs.filter_zero_coverage_lines && extension == '.pileup')
                          return " | awk '{if($4 != 0) print $0}'  > " + new_filename + extension
                      else
                          return " > " + new_filename + extension
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerImageId: 2fb927277493
                dockerPull: 'images.sbgenomics.com/milana_kaljevic/samtools:1.6'
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - SAM/BAM-Processing
            'sbg:cmdPreview': >-
              /opt/samtools-1.6/samtools mpileup    input_bam_file1.bam 
              input_bam_file2.bam   | awk '{if($4 != 0) print $0}'  >
              input_bam_file1_chr20.pileup
            'sbg:image_url': null
            'sbg:license': The MIT License
            'sbg:links':
              - id: 'http://samtools.sourceforge.net/'
                label: Homepage
              - id: 'https://github.com/samtools/samtools'
                label: Source code
              - id: 'http://sourceforge.net/p/samtools/wiki/Home/'
                label: Wiki
              - id: 'http://sourceforge.net/projects/samtools/files/'
                label: Download
              - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                label: Publication
              - id: 'http://www.htslib.org/doc/samtools.html'
                label: Documentation
            'sbg:toolAuthor': 'Heng Li, Sanger Institute'
            'sbg:toolkit': SAMtools
            'sbg:toolkitVersion': '1.6'
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510090
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SAMtools-Mpileup-1-6-fix/1
              - 'sbg:revision': 1
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510842
                'sbg:revisionNotes': 'Remove defaults, remove ''fix'' from tool name'
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/samtools-mpileup-1-6/1
            'sbg:revision': 1
            'sbg:revisionNotes': 'Remove defaults, remove ''fix'' from tool name'
            'sbg:modifiedOn': 1560510842
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510090
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 1
            'sbg:publisher': sbg
            'sbg:content_hash': a61bb5f8ff32335736d4ecaacdba6ee445f8ba77d77ff9e323ac5b52b3c106b98
          label: SAMtools Mpileup 1.6
          scatter:
            - region_pileup_generation
          'sbg:x': -137
          'sbg:y': 88
        - id: sbg_samtools_merge_mpileup
          in:
            - id: contig_order_names
              source:
                - sbg_prepare_intervals/names
            - id: file_list
              source:
                - samtools_mpileup_1_6/output_pileup_vcf_or_bcf_file
            - id: output_state
              default: Merge
          out:
            - id: output_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/sbg-samtools-merge-mpileup/1
            baseCommand: []
            inputs:
              - 'sbg:includeInPorts': true
                id: contig_order_names
                type: 'string[]'
                label: merging_order
                doc: >-
                  Ordered names of contigs that are used for merging of file
                  list loaded via file_list input.
              - id: file_list
                type: 'File[]'
                label: Input file list
                doc: Input files for merging.
                'sbg:fileTypes': 'PILEUP, VCF, BCF'
              - id: output_state
                type:
                  type: enum
                  symbols:
                    - Merge
                    - Pass nonempty files
                    - Pass all files
                  name: output_state
                label: Output state
                doc: >-
                  Mode to tell how to merge mpileup files or to pass them
                  through.
            outputs:
              - id: output_file
                doc: 'File obtained by merging files form #file_list.'
                label: output_file
                type: File?
                outputBinding:
                  glob: |-
                    ${
                        var new_files = [];
                        for (i = 0; i < inputs.file_list.length; i++) {
                            if (inputs.file_list[i] != null) {
                                new_files.push(inputs.file_list[i]);
                            }
                        }

                        if (inputs.output_state == "Merge") {

                            format = new_files[0].path.split("/").pop();
                            format = format.split('.').pop();

                            return "*.merged." + format.toLowerCase();
                        } else {
                            return "{*.pileup,*.vcf,*.bcf}";
                        }

                    }
                  outputEval: |-
                    ${

                        var inputs = {};

                        to_ret = "";

                        for (fn in inputs.file_list) {
                            if (inputs.file_list[fn].size !== 0) {

                                name = inputs.file_list[fn].path.split("/").pop();


                                if ("metadata" in inputs.file_list[fn]) {

                                    inputs[name] = inputs.file_list[fn].metadata;

                                }

                            }
                        }

                        content = [].concat(self);

                        if (inputs.output_state == "Merge") {

                            // set metadata from any input to this tool
                            name = content[0]["name"];

                            if (name) {
                                metadata = inputs[name]

                                if (metadata) {
                                    if ("ScatteredUsing" in metadata)
                                        delete metadata["ScatteredUsing"];

                                    for (var fp in content) {

                                        content[fp]["metadata"] = metadata;

                                    }
                                }
                            }

                        } else {
                            for (var fp in content) {

                                name = content[fp]["name"]

                                if (name in inputs) {
                                    metadata = inputs[name];

                                    content[fp]["metadata"] = metadata;
                                }

                            }

                            return content[0];
                        }

                        return content[0];
                    }
                'sbg:fileTypes': 'PILEUP, BCF, VCF'
            doc: >-
              Main function of this tool is to provide option to the user
              whether to merge list of files from SAMtools Mpileup tool or to
              pass list of files for each contig to the output. 


              If MERGE option is selected, regardless of the type, files will be
              merged following the contig order from the reference file index
              file or BED file. This order is provided using input
              contig_order_names.  So, this tool will receive list of pileup,
              vcf or BC files and list of files names of BED files used in
              creating input files. Pairing and finding the correct order will
              be done using metadata field  ScatteredUsing.


              ####Common issue:

              This tools is built specially for mpileup parallel and Varscan
              pipelines, so please don't use it outside these pipeline unless
              you are sure how to use it. Contact our support team for
              information.
            label: SBG SAMtools Merge Mpileup
            arguments:
              - position: 3
                prefix: ''
                shellQuote: false
                valueFrom: |-
                  ${
                      var new_files = [];
                      for (i = 0; i < inputs.file_list.length; i++) {
                          if (inputs.file_list[i] != null) {
                              new_files.push(inputs.file_list[i]);
                          }
                      }

                      if (inputs.output_state == "Merge") {
                          to_ret = "";
                          for (f in inputs.contig_order_names) {

                              for (fn in new_files) {
                                  if (new_files[fn].metadata) {

                                      if (inputs.contig_order_names[f] == new_files[fn].metadata.ScatteredUsing) {
                                          if (to_ret == "")
                                              to_ret = new_files[fn].path;
                                          else
                                              to_ret = to_ret + "," + new_files[fn].path;

                                      }
                                  }
                              }
                          }

                          format = new_files[0].path.split("/").pop();
                          format = format.split('.').pop();

                          if (to_ret == "")
                              return;
                          else
                              return "python merge.py --output-format " + format.toLowerCase() + " --file-list " + to_ret;

                      } else if (inputs.output_state == "Pass nonempty files") {
                          to_ret = "";

                          for (fn in new_files) {
                              if (new_files[fn].size == 0) {

                                  if (to_ret == "")
                                      to_ret = new_files[fn].path;
                                  else
                                      to_ret = to_ret + "," + new_files[fn].path;
                              }
                          }

                          if (to_ret == "")
                              return;
                          else
                              return "python delete_empty.py --file-list " + to_ret;
                      }
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'images.sbgenomics.com/marouf/sbg-pileup-merge-out:0.1'
              - class: InitialWorkDirRequirement
                listing:
                  - entryname: merge.py
                    entry: >-
                      """

                      Usage: merge.py     --file-list FILE --output-format
                      FRMT  


                      Description:
                        Main function of this tool is to provide option to the user weather to merge list of files from SAMtools Mpileup
                        tool or to pass list of files for each contig to the output. If --Merge option is selected, regardless of the type,
                        files will be merged following the contig order from the BAM header.

                      Options:

                          --help                              This message.

                          --version                           Tool version.
                          
                          --output-format FRMT                Output_format.

                          --file-list FILE                    List of contigs in BCF or PILEUP format produced by SAMtools Mpileup tool.

                      """


                      import docopt

                      from subprocess import Popen


                      args = docopt.docopt(__doc__, version="1.0") 


                      BCFTOOLS_ROOT = "bcftools"


                      sorted_file_list = list()

                      file_list = args['--file-list'].split(',')


                      output_format =  args['--output-format']



                      for fp in file_list:
                          sorted_file_list.append(fp)



                      if output_format =="pileup":
                          with open(sorted_file_list[0].split('/')[-1].replace('.pileup', '.merged.pileup'), 'w')  as rr:
                              cmd = ['cat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                      elif output_format =="bcf":
                          with open(sorted_file_list[0].split('/')[-1].replace('.bcf', '.merged.bcf'), 'w')  as rr:
                              cmd = [BCFTOOLS_ROOT, 'concat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                      elif output_format =="vcf":
                          with open(sorted_file_list[0].split('/')[-1].replace('.vcf', '.merged.vcf'), 'w')  as rr:
                              cmd = [BCFTOOLS_ROOT, 'concat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                    writable: false
                  - entryname: delete_empty.py
                    entry: |-
                      """
                      Usage: delete_empty.py     --file-list FILE  

                      Description:
                        Main function of this tool is to remove empty pileup files.

                      Options:

                          --help                              This message.

                          --version                           Tool version.

                          --file-list FILE                    List of contigs in BCF or PILEUP format produced by SAMtools Mpileup tool.

                      """

                      import os
                      import docopt
                      from subprocess import Popen

                      args = docopt.docopt(__doc__, version="1.0") 


                      sorted_file_list = list()
                      file_list = args['--file-list'].split(',')

                      for fp in file_list:
                          os.remove(fp)
                    writable: false
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - Text-Processing
              - SAM/BAM-Processing
            'sbg:cmdPreview': >-
              python merge.py --output-format vcf --file-list
              1994060146_RNASeq_R.Aligned.toTranscriptome.out_11.vcf,/ovde/1994060146_RNASeq_R.Aligned.toTranscriptome.out_22.vcf
            'sbg:image_url': null
            'sbg:license': Apache License 2.0
            'sbg:toolAuthor': >-
              Mohamed Marouf, Seven Bridges Genomics,
              <mohamed.marouf@sbgenomics.com>
            'sbg:toolkit': SBGTools
            'sbg:toolkitVersion': '0.2'
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510115
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SBG-SAMtools-Merge-Mpileup/2
              - 'sbg:revision': 1
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560512063
                'sbg:revisionNotes': Remove unnecessary create file
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/sbg-samtools-merge-mpileup/1
            'sbg:revision': 1
            'sbg:revisionNotes': Remove unnecessary create file
            'sbg:modifiedOn': 1560512063
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510115
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 1
            'sbg:publisher': sbg
            'sbg:content_hash': a2ec217f3a73795085cf637b4eec9bd3c01046c4c642ae76a554dac977b1b3992
          label: SBG SAMtools Merge Mpileup
          'sbg:x': 147.2751007080078
          'sbg:y': 30.081085205078125
      requirements:
        - class: ScatterFeatureRequirement
      'sbg:projectName': BMS CNV v2 - Dev
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': vtomic
          'sbg:modifiedOn': 1560512198
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': vtomic
          'sbg:modifiedOn': 1560512973
          'sbg:revisionNotes': Initial version
        - 'sbg:revision': 2
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1560516941
          'sbg:revisionNotes': added PILEUP as output
        - 'sbg:revision': 3
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1561052241
          'sbg:revisionNotes': exposed reference
      'sbg:image_url': >-
        https://igor.sbgenomics.com/ns/brood/images/vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3.png
      'sbg:appVersion':
        - v1.0
      'sbg:id': vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3
      'sbg:revision': 3
      'sbg:revisionNotes': exposed reference
      'sbg:modifiedOn': 1561052241
      'sbg:modifiedBy': ana_popic
      'sbg:createdOn': 1560512198
      'sbg:createdBy': vtomic
      'sbg:project': vtomic/bms-cnv-v2-dev
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
      'sbg:contributors':
        - ana_popic
        - vtomic
      'sbg:latestRevision': 3
      'sbg:publisher': sbg
      'sbg:content_hash': ab6b0ada5489223c262d31fb31b2d5ac543917f1f2e338e34f89e1a2a72bfa081
    label: SAMtools Mpileup Parallel 1.6 - Tumor
    'sbg:x': -226
    'sbg:y': -513
  - id: samtools_mpileup_parallel_1_6_normal
    in:
      - id: fai_file
        source: filtered_reference_fai
      - id: input_bam_files
        source:
          - normal_bam
      - id: reference_fasta
        source: reference_fasta
    out:
      - id: output_file
    run:
      class: Workflow
      cwlVersion: v1.0
      id: vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3
      label: SAMtools Mpileup Parallel 1.6
      $namespaces:
        sbg: 'https://sevenbridges.com'
      inputs:
        - id: fai_file
          'sbg:fileTypes': FAI
          type: File?
          label: Reference FAI
          'sbg:x': -509
          'sbg:y': -119
        - id: input_bam_files
          'sbg:fileTypes': BAM
          type: 'File[]'
          label: BAM files
          'sbg:x': -512
          'sbg:y': 150
        - id: reference_fasta
          'sbg:fileTypes': 'FASTA, FA, GZ'
          type: File?
          label: Reference FASTA
          'sbg:x': -508
          'sbg:y': 18
      outputs:
        - id: output_file
          outputSource:
            - sbg_samtools_merge_mpileup/output_file
          'sbg:fileTypes': 'PILEUP, BCF, VCF'
          type: File?
          label: Output file
          'sbg:x': 373.5552978515625
          'sbg:y': 27.004995346069336
      steps:
        - id: sbg_prepare_intervals
          in:
            - id: fai_file
              source: fai_file
            - id: format
              default: chr start end
            - id: split_mode
              default: File per chr with alt contig in a single file
          out:
            - id: intervals
            - id: names
            - id: str_arr
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/sbg-prepare-intervals/0
            baseCommand:
              - python
              - sbg_prepare_intervals.py
            inputs:
              - 'sbg:category': File Input
                id: bed_file
                type: File?
                inputBinding:
                  position: 1
                  prefix: '--bed'
                  shellQuote: false
                label: Input BED file
                doc: >-
                  Input BED file containing intervals. Required for modes 3 and
                  4.
                'sbg:fileTypes': BED
              - 'sbg:category': File Input
                id: fai_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '--fai'
                  shellQuote: false
                label: Input FAI file
                doc: >-
                  FAI file is converted to BED format if BED file is not
                  provided.
                'sbg:fileTypes': FAI
              - 'sbg:category': Input
                id: format
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - chr start end
                      - 'chr:start-end'
                    name: format
                label: Interval format
                doc: Format of the intervals in the generated files.
              - default: 0
                'sbg:category': Input
                id: split_mode
                type:
                  type: enum
                  symbols:
                    - File per interval
                    - File per chr with alt contig in a single file
                    - Output original BED
                    - File per interval with alt contig in a single file
                  name: split_mode
                inputBinding:
                  position: 3
                  prefix: '--mode'
                  shellQuote: false
                  valueFrom: |-
                    ${
                        if (self == 0) {
                            self = null;
                            inputs.split_mode = null
                        };


                        mode = inputs.split_mode
                        switch (mode) {
                            case "File per interval":
                                return 1
                            case "File per chr with alt contig in a single file":
                                return 2
                            case "Output original BED":
                                return 3
                            case "File per interval with alt contig in a single file":
                                return 4
                        }
                        return 3
                    }
                label: Split mode
                doc: >-
                  Depending on selected Split Mode value, output files are
                  generated in accordance with description below:  1. File per
                  interval - The tool creates one interval file per line of the
                  input BED(FAI) file. Each interval file contains a single line
                  (one of the lines of BED(FAI) input file).  2. File per chr
                  with alt contig in a single file - For each contig(chromosome)
                  a single file is created containing all the intervals
                  corresponding to it . All the intervals (lines) other than
                  (chr1, chr2 ... chrY or 1, 2 ... Y) are saved as
                  ("others.bed").  3. Output original BED - BED file is required
                  for execution of this mode. If mode 3 is applied input is
                  passed to the output.  4. File per interval with alt contig in
                  a single file - For each chromosome a single file is created
                  for each interval. All the intervals (lines) other than (chr1,
                  chr2 ... chrY or 1, 2 ... Y) are saved as ("others.bed").
                  NOTE: Do not use option 1 (File per interval) with exome BED
                  or a BED with a lot of GL contigs, as it will create a large
                  number of files.
            outputs:
              - id: intervals
                doc: Array of BED files genereted as per selected Split Mode.
                label: Intervals
                type: 'File[]?'
                outputBinding:
                  glob: Intervals/*.bed
                  outputEval: |-
                    ${

                        for (var i = 0; i < self.length; i++) {
                            var out_metadata = {
                                'sbg_scatter': 'true'
                            };
                            self[i] = setMetadata(self[i], out_metadata)
                        };

                        return self

                    }
                'sbg:fileTypes': BED
              - id: names
                doc: File containing the names of created files.
                label: Output file names
                type: string?
                outputBinding:
                  loadContents: true
                  glob: Intervals/names.txt
                  outputEval: |-
                    ${
                        content = self[0].contents.replace(/\0/g, '')
                        content = content.replace('[', '')
                        content = content.replace(']', '')
                        content = content.replace(/\'/g, "")
                        content = content.replace(/\s/g, '')
                        content_arr = content.split(",")

                        return content_arr


                    }
              - id: str_arr
                doc: Outputs BED content as strings.
                label: String output
                type: 'string[]?'
                outputBinding:
                  loadContents: true
                  glob: |-
                    ${
                        if (inputs.bed_file) {
                            glob = inputs.bed_file.path
                            glob = glob.split('/').slice(-1)[0]
                        } else if (inputs.fai_file) {
                            glob = inputs.fai_file.path
                            glob = glob.split('/').slice(-1)[0].split('.').slice(0, -1).join('.') + '.bed'
                        }

                        return glob
                    }
                  outputEval: |-
                    ${
                        rows = self[0].contents
                        if (rows[rows.length - 1] == '\n') {
                            rows = rows.split(/\r?\n/).slice(0, -1);
                        } else {
                            rows = rows.split(/\r?\n/);
                        }
                        out_list = []
                        for (i = 0; i < rows.length; i++) {
                            row = rows[i];
                            chromosome = row.split("\t")[0];
                            start = row.split("\t")[1];
                            end = row.split("\t")[2];
                            if (start) {
                                interval = chromosome.concat(":", start, "-", end);
                            } else {
                                interval = chromosome
                            }
                            out_list.push(interval);
                        }
                        return out_list;

                    }
            doc: >-
              Depending on selected Split Mode value, output files are generated
              in accordance with description below:


              1. File per interval - The tool creates one interval file per line
              of the input BED(FAI) file.

              Each interval file contains a single line (one of the lines of
              BED(FAI) input file).


              2. File per chr with alt contig in a single file - For each
              contig(chromosome) a single file

              is created containing all the intervals corresponding to it .

              All the intervals (lines) other than (chr1, chr2 ... chrY or 1, 2
              ... Y) are saved as

              ("others.bed").


              3. Output original BED - BED file is required for execution of
              this mode. If mode 3 is applied input is passed to the output.


              4. File per interval with alt contig in a single file - For each
              chromosome a single file is created for each interval.

              All the intervals (lines) other than (chr1, chr2 ... chrY or 1, 2
              ... Y) are saved as

              ("others.bed").


              ##### Common issues: 

              Do not use option 1 (File per interval) with exome BED or a BED
              with a lot of GL contigs, as it will create a large number of
              files.
            label: SBG Prepare Intervals
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: |-
                  ${
                      if (inputs.format)
                          return "--format " + "\"" + inputs.format + "\""
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'images.sbgenomics.com/bogdang/sbg_prepare_intervals:1.0'
              - class: InitialWorkDirRequirement
                listing:
                  - entryname: sbg_prepare_intervals.py
                    entry: >-
                      """

                      Usage:
                          sbg_prepare_intervals.py [options] [--fastq FILE --bed FILE --mode INT --format STR --others STR]

                      Description:
                          Purpose of this tool is to split BED file into files based on the selected mode.
                          If bed file is not provided fai(fasta index) file is converted to bed.

                      Options:

                          -h, --help            Show this message.

                          -v, -V, --version     Tool version.

                          -b, -B, --bed FILE    Path to input bed file.

                          --fai FILE            Path to input fai file.

                          --format STR          Output file format.

                          --mode INT            Select input mode.

                      """



                      import os

                      import sys

                      import glob

                      import shutil

                      from docopt import docopt


                      default_extension = '.bed'  # for output files



                      def create_file(contents, contig_name,
                      extension=default_extension):
                          """function for creating a file for all intervals in a contig"""

                          new_file = open("Intervals/" + contig_name + extension, "w")
                          new_file.write(contents)
                          new_file.close()


                      def add_to_file(line, name, extension=default_extension):
                          """function for adding a line to a file"""

                          new_file = open("Intervals/" + name + extension, "a")
                          if lformat == formats[1]:
                              sep = line.split("\t")
                              line = sep[0] + ":" + sep[1] + "-" + sep[2]
                          new_file.write(line)
                          new_file.close()


                      def fai2bed(fai):
                          """function to create a bed file from fai file"""

                          region_thr = 10000000  # threshold used to determine starting point accounting for telomeres in chromosomes
                          basename = fai[0:fai.rfind(".")]
                          with open(fai, "r") as ins:
                              new_array = []
                              for line in ins:
                                  len_reg = int(line.split()[1])
                                  cutoff = 0 if (
                                  len_reg < region_thr) else 0  # sd\\telomeres or start with 1
                                  new_line = line.split()[0] + '\t' + str(cutoff) + '\t' + str(
                                      len_reg + cutoff)
                                  new_array.append(new_line)
                          new_file = open(basename + ".bed", "w")
                          new_file.write("\n".join(new_array))
                          return basename + ".bed"


                      def chr_intervals(no_of_chrms=23):
                          """returns all possible designations for chromosome intervals"""

                          chrms = []
                          for i in range(1, no_of_chrms):
                              chrms.append("chr" + str(i))
                              chrms.append(str(i))
                          chrms.extend(["x", "y", "chrx", "chry"])
                          return chrms


                      def mode_1(orig_file):
                          """mode 1: every line is a new file"""

                          with open(orig_file, "r") as ins:
                              prev = ""
                              counter = 0
                              names = []
                              for line in ins:
                                  if is_header(line):
                                      continue
                                  if line.split()[0] == prev:
                                      counter += 1
                                  else:
                                      counter = 0
                                  suffix = "" if (counter == 0) else "_" + str(counter)
                                  create_file(line, line.split()[0] + suffix)
                                  names.append(line.split()[0] + suffix)
                                  prev = line.split()[0]

                              create_file(str(names), "names", extension=".txt")


                      def mode_2(orig_file, others_name):
                          """mode 2: separate file is created for each chromosome, and one file is created for other intervals"""

                          chrms = chr_intervals()
                          names = []

                          with open(orig_file, 'r') as ins:
                              for line in ins:
                                  if is_header(line):
                                      continue
                                  name = line.split()[0]
                                  if name.lower() in chrms:
                                      name = name
                                  else:
                                      name = others_name
                                  try:
                                      add_to_file(line, name)
                                      if not name in names:
                                          names.append(name)
                                  except:
                                      raise Exception(
                                          "Couldn't create or write in the file in mode 2")

                              create_file(str(names), "names", extension=".txt")


                      def mode_3(orig_file, extension=default_extension):
                          """mode 3: input file is staged to output"""

                          orig_name = orig_file.split("/")[len(orig_file.split("/")) - 1]
                          output_file = r"./Intervals/" + orig_name[
                                                          0:orig_name.rfind('.')] + extension

                          shutil.copyfile(orig_file, output_file)

                          names = [orig_name[0:orig_name.rfind('.')]]
                          create_file(str(names), "names", extension=".txt")


                      def mode_4(orig_file, others_name):
                          """mode 4: every interval in chromosomes is in a separate file. Other intervals are in a single file"""

                          chrms = chr_intervals()
                          names = []

                          with open(orig_file, "r") as ins:
                              counter = {}
                              for line in ins:
                                  if line.startswith('@'):
                                      continue
                                  name = line.split()[0].lower()
                                  if name in chrms:
                                      if name in counter:
                                          counter[name] += 1
                                      else:
                                          counter[name] = 0
                                      suffix = "" if (counter[name] == 0) else "_" + str(counter[name])
                                      create_file(line, name + suffix)
                                      names.append(name + suffix)
                                      prev = name
                                  else:
                                      name = others_name
                                      if not name in names:
                                          names.append(name)
                                      try:
                                          add_to_file(line, name)
                                      except:
                                          raise Exception(
                                              "Couldn't create or write in the file in mode 4")

                          create_file(str(names), "names", extension=".txt")


                      def prepare_intervals():
                          # reading input files and split mode from command line
                          args = docopt(__doc__, version='1.0')

                          bed_file = args['--bed']
                          fai_file = args['--fai']
                          split_mode = int(args['--mode'])

                          # define file name for non-chromosomal contigs
                          others_name = 'others'

                          global formats, lformat
                          formats = ["chr start end", "chr:start-end"]
                          lformat = args['--format']
                          if lformat == None:
                              lformat = formats[0]
                          if not lformat in formats:
                              raise Exception('Unsuported interval format')

                          if not os.path.exists(r"./Intervals"):
                              os.mkdir(r"./Intervals")
                          else:
                              files = glob.glob(r"./Intervals/*")
                              for f in files:
                                  os.remove(f)

                          # create variable input_file taking bed_file as priority
                          if bed_file:
                              input_file = bed_file
                          elif fai_file:
                              input_file = fai2bed(fai_file)
                          else:
                              raise Exception('No input files are provided')

                          # calling adequate split mode function
                          if split_mode == 1:
                              mode_1(input_file)
                          elif split_mode == 2:
                              mode_2(input_file, others_name)
                          elif split_mode == 3:
                              if bed_file:
                                  mode_3(input_file)
                              else:
                                  raise Exception('Bed file is required for mode 3')
                          elif split_mode == 4:
                              mode_4(input_file, others_name)
                          else:
                              raise Exception('Split mode value is not set')


                      def is_header(line):
                          x = line.split('\t')
                          try:
                              int(x[1])
                              int(x[2])
                              header = False
                          except:
                              sys.stderr.write('Line is skipped: {}'.format(line))
                              header = True
                          return header


                      if __name__ == '__main__':
                          prepare_intervals()
                    writable: false
                  - $(inputs.fai_file)
                  - $(inputs.bed_file)
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - Converters
            'sbg:cmdPreview': python sbg_prepare_intervals.py  --format "chr start end" --mode 2
            'sbg:image_url': null
            'sbg:license': Apache License 2.0
            'sbg:toolAuthor': Seven Bridges Genomics
            'sbg:toolkit': SBGTools
            'sbg:toolkitVersion': '1.0'
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/sbg-prepare-intervals/0
            'sbg:revision': 0
            'sbg:revisionNotes': >-
              Copy of
              vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
            'sbg:modifiedOn': 1560510063
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510063
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 0
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510063
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
            'sbg:publisher': sbg
            'sbg:content_hash': aa51632b14cdb5f6e798ed50dbb067444bf78dfa99bbdb384445c93e9f8782262
            'sbg:copyOf': >-
              vojislav_varjacic/vojislav-varjacics-demo-project/SBG-Prepare-Intervals/0
          label: SBG Prepare Intervals
          'sbg:x': -294
          'sbg:y': -63
        - id: samtools_mpileup_1_6
          in:
            - id: input_bam_files
              source:
                - input_bam_files
            - id: output_format
              default: PILEUP
            - id: reference_fasta_file
              source: reference_fasta
            - id: region_pileup_generation
              source: sbg_prepare_intervals/names
          out:
            - id: output_pileup_vcf_or_bcf_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/samtools-mpileup-1-6/1
            baseCommand: []
            inputs:
              - 'sbg:category': Input options
                id: assume_the_quality_Illumina_encoding
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-6'
                  shellQuote: false
                label: Assume the quality is in the Illumina-1.3+ encoding
                doc: Assume the quality is in the Illumina-1.3+ encoding.
              - 'sbg:category': File input
                id: chr_pos_list_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-l'
                  shellQuote: false
                label: List of positions (chr pos) or regions in BED file
                doc: >-
                  BED or position list file containing a list of regions or
                  sites where pileup or BCF should be generated.
                'sbg:fileTypes': BED
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: comma_separated_list_of_platforms_for_indels
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-P'
                  shellQuote: false
                label: Comma separated list of platforms for indels
                doc: >-
                  Comma dilimited list of platforms (determined by @RG-PL) from
                  which indel candidates are obtained. It is recommended to
                  collect indel candidates from sequencing technologies that
                  have low indel error rate such as ILLUMINA.
              - 'sbg:category': Input options
                id: count_anomalous_read_pairs
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-A'
                  shellQuote: false
                label: Count anomalous read pairs
                doc: Do not skip anomalous read pairs in variant calling.
              - 'sbg:category': Input options
                id: disable_baq_computation
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-B'
                  shellQuote: false
                label: Disable BAQ computation
                doc: >-
                  Disable probabilistic realignment for the computation of base
                  alignment quality (BAQ). BAQ is the Phred-scaled probability
                  of a read base being misaligned. Applying this option greatly
                  helps to reduce false SNPs caused by misalignments.
              - 'sbg:category': File input
                id: exclude_read_groups_list_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-G'
                  shellQuote: false
                label: Exclude read groups listed in file
                doc: Exclude read groups listed in file.
                'sbg:fileTypes': TXT
              - 'sbg:category': Input options
                id: filter_flags
                type: string?
                inputBinding:
                  position: 2
                  prefix: '--ff'
                  shellQuote: false
                label: Filter flags
                doc: >-
                  Filter flags: skip reads with mask bits set
                  [UNMAP,SECONDARY,QCFAIL,DUP].
              - 'sbg:category': output configuration
                id: filter_zero_coverage_lines
                type: boolean?
                label: Filter zero coverage lines.
                doc: >-
                  Filter zero coverage lines from pileup files. This option is
                  valid only when output pileup files.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: gap_extension_seq_error_probability
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-e'
                  shellQuote: false
                label: Phred-scaled gap extension seq error probability
                doc: >-
                  Phred-scaled gap extension sequencing error probability.
                  Reducing INT leads to longer indels.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: gap_open_sequencing_error_probability
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-o'
                  shellQuote: false
                label: Phred-scaled gap open sequencing error probability
                doc: >-
                  Phred-scaled gap open sequencing error probability. Reducing
                  INT leads to more indel calls.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: homopolymer_err_coeficient
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-h'
                  shellQuote: false
                label: Coefficient for homopolymer errors
                doc: >-
                  Coefficient for modeling homopolymer errors. Given an l-long
                  homopolymer run, the sequencing error of an indel of size s is
                  modeled as INT*s/l.
              - 'sbg:category': configuration
                id: ignore_overlaps
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-x'
                  shellQuote: false
                label: Disable read-pair overlap detection
                doc: Disable read-pair overlap detection.
              - 'sbg:category': Input options
                id: ignore_rg_tags
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-R'
                  shellQuote: false
                label: Ignore RG tags
                doc: Ignore rg tags.
              - 'sbg:category': Input options
                schema:
                  - File
                id: input_bam_files
                type: 'File[]'
                inputBinding:
                  position: 3
                  shellQuote: false
                  separator: ' '
                  streamable: false
                label: List of input BAM files
                doc: List of input BAM files.
                'sbg:fileTypes': BAM
                secondaryFiles:
                  - ^.bai
              - 'sbg:category': Input options
                id: mapq_threshold
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-q'
                  shellQuote: false
                  separator: ' '
                label: Skip alignments with mapQ  smaller than
                doc: Minimum mapping quality for an alignment to be used.
              - 'sbg:category': configuration
                schema:
                  - 'null'
                  - int
                id: max_idepth
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-L'
                  shellQuote: false
                label: Skip INDEL calling if the average per-sample depth is above
                doc: Skip INDEL calling if the average per-sample depth is above.
              - 'sbg:category': Input options
                id: max_per_bam_depth
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-d'
                  shellQuote: false
                label: Max per-BAM depth
                doc: 'At a position, read maximally INT reads per input BAM.'
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: min_gapped_reads_for_indel
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-m'
                  shellQuote: false
                label: Minimum gapped reads for indel candidates
                doc: Minimum gapped reads for indel candidates.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: minimum_fraction_of_gapped_reads
                type: float?
                inputBinding:
                  position: 2
                  prefix: '-F'
                  shellQuote: false
                label: Minimum fraction of gapped reads for candidates
                doc: Minimum fraction of gapped reads for candidates.
              - 'sbg:category': configuration
                id: more_info_to_output
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-t'
                  shellQuote: false
                label: Comma-separated list of FORMAT and INFO tags to output
                doc: >-
                  Comma-separated list of FORMAT and INFO tags
                  (DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []") to output.
              - 'sbg:category': SNP/INDEL genotype likelihoods options
                id: no_indel_calling
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-I'
                  shellQuote: false
                label: Do not perform indel calling
                doc: Do not perform indel calling.
              - 'sbg:category': Output options
                id: output_base_positions_on_reads
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-O'
                  shellQuote: false
                label: Output base positions on reads
                doc: Output base positions on reads (disabled by -g/-u).
              - 'sbg:category': configuration
                id: output_format
                type:
                  type: enum
                  symbols:
                    - PILEUP
                    - BCF
                    - VCF
                  name: output_format
                label: Output file format
                doc: Output file format.
              - 'sbg:category': Output options
                id: output_mapping_quality
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-s'
                  shellQuote: false
                label: Output mapping quality
                doc: Output mapping quality (disabled by -g/-u).
              - 'sbg:category': Output options
                id: output_uncompressed_bcf_or_vcf
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-u'
                  shellQuote: false
                label: Generate uncompress BCF/VCF output
                doc: >-
                  Similar to -g except that the output is uncompressed BCF,
                  which is preferred for piping.
              - 'sbg:category': Input options
                id: parameter_for_adjusting_mapq
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-C'
                  shellQuote: false
                label: Parameter for adjusting mapQ
                doc: >-
                  Coefficient for downgrading mapping quality for reads
                  containing excessive mismatches. Given a read with a
                  phred-scaled probability q of being generated from the mapped
                  position, the new mapping quality is about
                  sqrt((INT-q)/INT)*INT. A zero value disables this
                  functionality; if enabled, the recommended value for BWA is
                  50.
              - 'sbg:category': configuration
                id: per_sample_mF
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-p'
                  shellQuote: false
                label: Apply -m and -F thresholds per sample
                doc: >-
                  Apply -m and -F thresholds per sample to increase sensitivity
                  of calling. By default both options are applied to reads
                  pooled from all samples.
              - 'sbg:category': Input options
                id: recalc_BAQ_fly
                type: boolean?
                inputBinding:
                  position: 2
                  prefix: '-E'
                  shellQuote: false
                  separator: ' '
                label: Recalculate BAQ on the fly
                doc: 'Recalculate BAQ on the fly, ignore existing BQ tags.'
              - 'sbg:category': File input
                schema:
                  - 'null'
                  - File
                id: reference_fasta_file
                type: File?
                inputBinding:
                  position: 2
                  prefix: '-f'
                  shellQuote: false
                label: Reference FASTA file
                doc: Faidx indexed reference sequence (FASTA) file.
                'sbg:fileTypes': 'FASTA, FA, GZ'
                secondaryFiles:
                  - .fai
              - 'sbg:category': Input options
                'sbg:includeInPorts': true
                id: region_pileup_generation
                type: string?
                inputBinding:
                  position: 2
                  prefix: '-r'
                  shellQuote: false
                label: Region in which pileup is generated
                doc: Region in which pileup is generated.
              - 'sbg:category': Input options
                id: required_flags
                type: int?
                inputBinding:
                  position: 2
                  prefix: '--rf'
                  shellQuote: false
                label: Required flags
                doc: 'Required flags: skip reads with mask bits unset.'
              - 'sbg:category': Input options
                id: skip_bases_with_baq_smaller_than_defined
                type: int?
                inputBinding:
                  position: 2
                  prefix: '-Q'
                  shellQuote: false
                  separator: ' '
                label: Skip bases with baseQ/BAQ smaller than
                doc: Minimum base quality for a base to be considered.
            outputs:
              - id: output_pileup_vcf_or_bcf_file
                doc: 'Output PILEUP, VCF, or BCF file.'
                label: 'Output PILEUP, VCF, or BCF file'
                type: File?
                outputBinding:
                  glob: '{*.pileup,*.bcf,*.vcf}'
                  outputEval: |-
                    ${
                        self = inheritMetadata(self, inputs.input_bam_files);


                        var add_metadata_key_ScatteredUsing = function(self, inputs) {
                            if (inputs.chr_pos_list_file) {
                                len = inputs.chr_pos_list_file.path.split('/').length

                                return inputs.chr_pos_list_file.path.split('/')[len - 1].split('.')[0]
                            } else if (inputs.region_pileup_generation) {
                                return inputs.region_pileup_generation
                            } else {
                                return
                            }
                        }
                        for (var i = 0; i < self.length; i++) {
                            var out_metadata = {
                                'ScatteredUsing': add_metadata_key_ScatteredUsing(self[i], inputs)
                            };
                            self[i] = setMetadata(self[i], out_metadata)
                        };

                        return self

                    }
                'sbg:fileTypes': 'PILEUP, BCF, VCF'
            doc: >-
              SAMtools Mpileup generates BCF or PILEUP for one or multiple BAM
              files. Alignment records are grouped by sample identifiers in @RG
              header lines. If sample identifiers are absent, each input file is
              regarded as one sample.


              In the pileup format (without -uor-g), each line represents a
              genomic position consisting of chromosome name, coordinate,
              reference base, read bases, read qualities, and alignment mapping
              qualities. Information on match, mismatch, indel, strand, mapping
              quality, and the start and end of a read are all encoded at the
              read base column. In this column, a dot stands for a match to the
              reference base on the forward strand, a comma for a match on the
              reverse strand, a '>' or '<' for a reference skip, 'ACGTN' for a
              mismatch on the forward strand, and 'acgtn' for a mismatch on the
              reverse strand. A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there
              is an insertion between this reference position and the next
              reference position. The length of the insertion is given by the
              integer in the pattern followed by the inserted sequence.
              Similarly, a pattern '-[0-9]+[ACGTNacgtn]+' represents a deletion
              from the reference. The deleted bases will be presented as `*' in
              the following lines. Also at the read base column, a symbol '^',
              marks the start of a read. The ASCII of the character following
              '^' minus 33 gives the mapping quality. The symbol '$' marks the
              end of a read segment.


              #### Common Issue:

              - Please use the public pileup parallel pipeline for large input
              files.

              - Please sort BAM files by coordinates before using mpileup.
            label: SAMtools Mpileup 1.6
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: /opt/samtools-1.6/samtools
              - position: 1
                shellQuote: false
                valueFrom: mpileup
              - position: 2
                shellQuote: false
                valueFrom: |-
                  ${

                      if (inputs.output_format) {
                          if (inputs.output_format == "PILEUP")
                              return
                          else if (inputs.output_format == "VCF")
                              return "-v"
                          else if (inputs.output_format == "BCF")
                              return "-g"
                          else
                              return


                      }



                  }
              - position: 1002
                shellQuote: false
                valueFrom: |-
                  ${

                      // Check if input files are delivered in an array or not
                      if (inputs.input_bam_files.constructor == Array) {
                          // Select the first entry of the array
                          filepath = "/" + inputs.input_bam_files[0].path
                      } else

                      {
                          // Use the file name directly
                          filepath = "/" + inputs.input_bam_files.path
                      }


                      filename = filepath.split("/").pop()
                      new_filename = filename.substr(0, filename.lastIndexOf("."))


                      if (inputs.chr_pos_list_file)
                          new_filename = new_filename + '_' + inputs.chr_pos_list_file.path.split('/')[inputs.chr_pos_list_file.path.split('/').length - 1].split('.')[0]

                      if (inputs.region_pileup_generation)
                          new_filename = new_filename + '_' + inputs.region_pileup_generation

                      extension = '.pileup'

                      if ((inputs.output_format && inputs.output_format === "BCF"))
                          extension = '.bcf'

                      else if ((inputs.output_format && inputs.output_format === "VCF"))
                          extension = '.vcf'

                      if (inputs.filter_zero_coverage_lines && extension == '.pileup')
                          return " | awk '{if($4 != 0) print $0}'  > " + new_filename + extension
                      else
                          return " > " + new_filename + extension
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerImageId: 2fb927277493
                dockerPull: 'images.sbgenomics.com/milana_kaljevic/samtools:1.6'
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - SAM/BAM-Processing
            'sbg:cmdPreview': >-
              /opt/samtools-1.6/samtools mpileup    input_bam_file1.bam 
              input_bam_file2.bam   | awk '{if($4 != 0) print $0}'  >
              input_bam_file1_chr20.pileup
            'sbg:image_url': null
            'sbg:license': The MIT License
            'sbg:links':
              - id: 'http://samtools.sourceforge.net/'
                label: Homepage
              - id: 'https://github.com/samtools/samtools'
                label: Source code
              - id: 'http://sourceforge.net/p/samtools/wiki/Home/'
                label: Wiki
              - id: 'http://sourceforge.net/projects/samtools/files/'
                label: Download
              - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                label: Publication
              - id: 'http://www.htslib.org/doc/samtools.html'
                label: Documentation
            'sbg:toolAuthor': 'Heng Li, Sanger Institute'
            'sbg:toolkit': SAMtools
            'sbg:toolkitVersion': '1.6'
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510090
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SAMtools-Mpileup-1-6-fix/1
              - 'sbg:revision': 1
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510842
                'sbg:revisionNotes': 'Remove defaults, remove ''fix'' from tool name'
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/samtools-mpileup-1-6/1
            'sbg:revision': 1
            'sbg:revisionNotes': 'Remove defaults, remove ''fix'' from tool name'
            'sbg:modifiedOn': 1560510842
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510090
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 1
            'sbg:publisher': sbg
            'sbg:content_hash': a61bb5f8ff32335736d4ecaacdba6ee445f8ba77d77ff9e323ac5b52b3c106b98
          label: SAMtools Mpileup 1.6
          scatter:
            - region_pileup_generation
          'sbg:x': -137
          'sbg:y': 88
        - id: sbg_samtools_merge_mpileup
          in:
            - id: contig_order_names
              source:
                - sbg_prepare_intervals/names
            - id: file_list
              source:
                - samtools_mpileup_1_6/output_pileup_vcf_or_bcf_file
            - id: output_state
              default: Merge
          out:
            - id: output_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: vtomic/bms-cnv-v2-dev/sbg-samtools-merge-mpileup/1
            baseCommand: []
            inputs:
              - 'sbg:includeInPorts': true
                id: contig_order_names
                type: 'string[]'
                label: merging_order
                doc: >-
                  Ordered names of contigs that are used for merging of file
                  list loaded via file_list input.
              - id: file_list
                type: 'File[]'
                label: Input file list
                doc: Input files for merging.
                'sbg:fileTypes': 'PILEUP, VCF, BCF'
              - id: output_state
                type:
                  type: enum
                  symbols:
                    - Merge
                    - Pass nonempty files
                    - Pass all files
                  name: output_state
                label: Output state
                doc: >-
                  Mode to tell how to merge mpileup files or to pass them
                  through.
            outputs:
              - id: output_file
                doc: 'File obtained by merging files form #file_list.'
                label: output_file
                type: File?
                outputBinding:
                  glob: |-
                    ${
                        var new_files = [];
                        for (i = 0; i < inputs.file_list.length; i++) {
                            if (inputs.file_list[i] != null) {
                                new_files.push(inputs.file_list[i]);
                            }
                        }

                        if (inputs.output_state == "Merge") {

                            format = new_files[0].path.split("/").pop();
                            format = format.split('.').pop();

                            return "*.merged." + format.toLowerCase();
                        } else {
                            return "{*.pileup,*.vcf,*.bcf}";
                        }

                    }
                  outputEval: |-
                    ${

                        var inputs = {};

                        to_ret = "";

                        for (fn in inputs.file_list) {
                            if (inputs.file_list[fn].size !== 0) {

                                name = inputs.file_list[fn].path.split("/").pop();


                                if ("metadata" in inputs.file_list[fn]) {

                                    inputs[name] = inputs.file_list[fn].metadata;

                                }

                            }
                        }

                        content = [].concat(self);

                        if (inputs.output_state == "Merge") {

                            // set metadata from any input to this tool
                            name = content[0]["name"];

                            if (name) {
                                metadata = inputs[name]

                                if (metadata) {
                                    if ("ScatteredUsing" in metadata)
                                        delete metadata["ScatteredUsing"];

                                    for (var fp in content) {

                                        content[fp]["metadata"] = metadata;

                                    }
                                }
                            }

                        } else {
                            for (var fp in content) {

                                name = content[fp]["name"]

                                if (name in inputs) {
                                    metadata = inputs[name];

                                    content[fp]["metadata"] = metadata;
                                }

                            }

                            return content[0];
                        }

                        return content[0];
                    }
                'sbg:fileTypes': 'PILEUP, BCF, VCF'
            doc: >-
              Main function of this tool is to provide option to the user
              whether to merge list of files from SAMtools Mpileup tool or to
              pass list of files for each contig to the output. 


              If MERGE option is selected, regardless of the type, files will be
              merged following the contig order from the reference file index
              file or BED file. This order is provided using input
              contig_order_names.  So, this tool will receive list of pileup,
              vcf or BC files and list of files names of BED files used in
              creating input files. Pairing and finding the correct order will
              be done using metadata field  ScatteredUsing.


              ####Common issue:

              This tools is built specially for mpileup parallel and Varscan
              pipelines, so please don't use it outside these pipeline unless
              you are sure how to use it. Contact our support team for
              information.
            label: SBG SAMtools Merge Mpileup
            arguments:
              - position: 3
                prefix: ''
                shellQuote: false
                valueFrom: |-
                  ${
                      var new_files = [];
                      for (i = 0; i < inputs.file_list.length; i++) {
                          if (inputs.file_list[i] != null) {
                              new_files.push(inputs.file_list[i]);
                          }
                      }

                      if (inputs.output_state == "Merge") {
                          to_ret = "";
                          for (f in inputs.contig_order_names) {

                              for (fn in new_files) {
                                  if (new_files[fn].metadata) {

                                      if (inputs.contig_order_names[f] == new_files[fn].metadata.ScatteredUsing) {
                                          if (to_ret == "")
                                              to_ret = new_files[fn].path;
                                          else
                                              to_ret = to_ret + "," + new_files[fn].path;

                                      }
                                  }
                              }
                          }

                          format = new_files[0].path.split("/").pop();
                          format = format.split('.').pop();

                          if (to_ret == "")
                              return;
                          else
                              return "python merge.py --output-format " + format.toLowerCase() + " --file-list " + to_ret;

                      } else if (inputs.output_state == "Pass nonempty files") {
                          to_ret = "";

                          for (fn in new_files) {
                              if (new_files[fn].size == 0) {

                                  if (to_ret == "")
                                      to_ret = new_files[fn].path;
                                  else
                                      to_ret = to_ret + "," + new_files[fn].path;
                              }
                          }

                          if (to_ret == "")
                              return;
                          else
                              return "python delete_empty.py --file-list " + to_ret;
                      }
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 1000
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'images.sbgenomics.com/marouf/sbg-pileup-merge-out:0.1'
              - class: InitialWorkDirRequirement
                listing:
                  - entryname: merge.py
                    entry: >-
                      """

                      Usage: merge.py     --file-list FILE --output-format
                      FRMT  


                      Description:
                        Main function of this tool is to provide option to the user weather to merge list of files from SAMtools Mpileup
                        tool or to pass list of files for each contig to the output. If --Merge option is selected, regardless of the type,
                        files will be merged following the contig order from the BAM header.

                      Options:

                          --help                              This message.

                          --version                           Tool version.
                          
                          --output-format FRMT                Output_format.

                          --file-list FILE                    List of contigs in BCF or PILEUP format produced by SAMtools Mpileup tool.

                      """


                      import docopt

                      from subprocess import Popen


                      args = docopt.docopt(__doc__, version="1.0") 


                      BCFTOOLS_ROOT = "bcftools"


                      sorted_file_list = list()

                      file_list = args['--file-list'].split(',')


                      output_format =  args['--output-format']



                      for fp in file_list:
                          sorted_file_list.append(fp)



                      if output_format =="pileup":
                          with open(sorted_file_list[0].split('/')[-1].replace('.pileup', '.merged.pileup'), 'w')  as rr:
                              cmd = ['cat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                      elif output_format =="bcf":
                          with open(sorted_file_list[0].split('/')[-1].replace('.bcf', '.merged.bcf'), 'w')  as rr:
                              cmd = [BCFTOOLS_ROOT, 'concat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                      elif output_format =="vcf":
                          with open(sorted_file_list[0].split('/')[-1].replace('.vcf', '.merged.vcf'), 'w')  as rr:
                              cmd = [BCFTOOLS_ROOT, 'concat'] + sorted_file_list
                              p = Popen(cmd, stdout=rr)
                              p.wait()
                    writable: false
                  - entryname: delete_empty.py
                    entry: |-
                      """
                      Usage: delete_empty.py     --file-list FILE  

                      Description:
                        Main function of this tool is to remove empty pileup files.

                      Options:

                          --help                              This message.

                          --version                           Tool version.

                          --file-list FILE                    List of contigs in BCF or PILEUP format produced by SAMtools Mpileup tool.

                      """

                      import os
                      import docopt
                      from subprocess import Popen

                      args = docopt.docopt(__doc__, version="1.0") 


                      sorted_file_list = list()
                      file_list = args['--file-list'].split(',')

                      for fp in file_list:
                          os.remove(fp)
                    writable: false
              - class: InlineJavascriptRequirement
                expressionLib:
                  - |-
                    var updateMetadata = function(file, key, value) {
                        file['metadata'][key] = value;
                        return file;
                    };


                    var setMetadata = function(file, metadata) {
                        if (!('metadata' in file)) {
                            file['metadata'] = {}
                        }
                        for (var key in metadata) {
                            file['metadata'][key] = metadata[key];
                        }
                        return file
                    };

                    var inheritMetadata = function(o1, o2) {
                        var commonMetadata = {};
                        if (!Array.isArray(o2)) {
                            o2 = [o2]
                        }
                        for (var i = 0; i < o2.length; i++) {
                            var example = o2[i]['metadata'];
                            for (var key in example) {
                                if (i == 0)
                                    commonMetadata[key] = example[key];
                                else {
                                    if (!(commonMetadata[key] == example[key])) {
                                        delete commonMetadata[key]
                                    }
                                }
                            }
                        }
                        if (!Array.isArray(o1)) {
                            o1 = setMetadata(o1, commonMetadata)
                        } else {
                            for (var i = 0; i < o1.length; i++) {
                                o1[i] = setMetadata(o1[i], commonMetadata)
                            }
                        }
                        return o1;
                    };

                    var toArray = function(file) {
                        return [].concat(file);
                    };

                    var groupBy = function(files, key) {
                        var groupedFiles = [];
                        var tempDict = {};
                        for (var i = 0; i < files.length; i++) {
                            var value = files[i]['metadata'][key];
                            if (value in tempDict)
                                tempDict[value].push(files[i]);
                            else tempDict[value] = [files[i]];
                        }
                        for (var key in tempDict) {
                            groupedFiles.push(tempDict[key]);
                        }
                        return groupedFiles;
                    };

                    var orderBy = function(files, key, order) {
                        var compareFunction = function(a, b) {
                            if (a['metadata'][key].constructor === Number) {
                                return a['metadata'][key] - b['metadata'][key];
                            } else {
                                var nameA = a['metadata'][key].toUpperCase();
                                var nameB = b['metadata'][key].toUpperCase();
                                if (nameA < nameB) {
                                    return -1;
                                }
                                if (nameA > nameB) {
                                    return 1;
                                }
                                return 0;
                            }
                        };

                        files = files.sort(compareFunction);
                        if (order == undefined || order == "asc")
                            return files;
                        else
                            return files.reverse();
                    };
            'sbg:categories':
              - Text-Processing
              - SAM/BAM-Processing
            'sbg:cmdPreview': >-
              python merge.py --output-format vcf --file-list
              1994060146_RNASeq_R.Aligned.toTranscriptome.out_11.vcf,/ovde/1994060146_RNASeq_R.Aligned.toTranscriptome.out_22.vcf
            'sbg:image_url': null
            'sbg:license': Apache License 2.0
            'sbg:toolAuthor': >-
              Mohamed Marouf, Seven Bridges Genomics,
              <mohamed.marouf@sbgenomics.com>
            'sbg:toolkit': SBGTools
            'sbg:toolkitVersion': '0.2'
            'sbg:projectName': BMS CNV v2 - Dev
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560510115
                'sbg:revisionNotes': >-
                  Copy of
                  vojislav_varjacic/vojislav-varjacics-demo-project/SBG-SAMtools-Merge-Mpileup/2
              - 'sbg:revision': 1
                'sbg:modifiedBy': vtomic
                'sbg:modifiedOn': 1560512063
                'sbg:revisionNotes': Remove unnecessary create file
            'sbg:appVersion':
              - v1.0
            'sbg:id': vtomic/bms-cnv-v2-dev/sbg-samtools-merge-mpileup/1
            'sbg:revision': 1
            'sbg:revisionNotes': Remove unnecessary create file
            'sbg:modifiedOn': 1560512063
            'sbg:modifiedBy': vtomic
            'sbg:createdOn': 1560510115
            'sbg:createdBy': vtomic
            'sbg:project': vtomic/bms-cnv-v2-dev
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - vtomic
            'sbg:latestRevision': 1
            'sbg:publisher': sbg
            'sbg:content_hash': a2ec217f3a73795085cf637b4eec9bd3c01046c4c642ae76a554dac977b1b3992
          label: SBG SAMtools Merge Mpileup
          'sbg:x': 147.2751007080078
          'sbg:y': 30.081085205078125
      requirements:
        - class: ScatterFeatureRequirement
      'sbg:projectName': BMS CNV v2 - Dev
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': vtomic
          'sbg:modifiedOn': 1560512198
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': vtomic
          'sbg:modifiedOn': 1560512973
          'sbg:revisionNotes': Initial version
        - 'sbg:revision': 2
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1560516941
          'sbg:revisionNotes': added PILEUP as output
        - 'sbg:revision': 3
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1561052241
          'sbg:revisionNotes': exposed reference
      'sbg:image_url': >-
        https://igor.sbgenomics.com/ns/brood/images/vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3.png
      'sbg:appVersion':
        - v1.0
      'sbg:id': vtomic/bms-cnv-v2-dev/samtools-mpileup-parallel-1-6/3
      'sbg:revision': 3
      'sbg:revisionNotes': exposed reference
      'sbg:modifiedOn': 1561052241
      'sbg:modifiedBy': ana_popic
      'sbg:createdOn': 1560512198
      'sbg:createdBy': vtomic
      'sbg:project': vtomic/bms-cnv-v2-dev
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
      'sbg:contributors':
        - ana_popic
        - vtomic
      'sbg:latestRevision': 3
      'sbg:publisher': sbg
      'sbg:content_hash': ab6b0ada5489223c262d31fb31b2d5ac543917f1f2e338e34f89e1a2a72bfa081
    label: SAMtools Mpileup Parallel 1.6 - Normal
    'sbg:x': -317
    'sbg:y': -183
  - id: control_freec_11_5
    in:
      - id: capture_regions
        source: capture_regions
      - id: chr_len
        source: chromosome_sizes
      - id: mate_file_control
        source: normal_bam
      - id: mate_file_sample
        source: tumor_bam
      - id: mate_orientation_control
        source: mate_orientation_control
      - id: mate_orientation_sample
        source: mate_orientation_sample
      - id: mini_pileup_control
        source: samtools_mpileup_parallel_1_6_normal/output_file
      - id: mini_pileup_sample
        source: samtools_mpileup_parallel_1_6_tumor/output_file
      - id: noisy_data
        source: noisy_data
      - id: ploidy
        source:
          - ploidy
      - id: reference
        source: reference_fasta
      - id: sex
        source: sex
      - id: snp_file
        source: snp_file
    out:
      - id: GC_profile
      - id: cnvs
      - id: cnvs_pvalue
      - id: config_script
      - id: control_cpn
      - id: control_pileup
      - id: info_txt
      - id: pngs
      - id: ratio
      - id: ratio_BedGraph
      - id: sample_BAF
      - id: sample_cpn
      - id: sample_pileup
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      id: vtomic/bms-cnv-v2-dev/control-freec-11-5/4
      baseCommand: []
      inputs:
        - 'sbg:category': General
          id: GC_content_profile
          type: File?
          label: GC content profile
          doc: GC-content profile for a given window-size.
          'sbg:fileTypes': CNP
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: bed_graph_output
          type: boolean?
          label: Bed Graph output
          doc: >-
            Set "BedGraphOutput=TRUE" if you want an additional output in
            BedGraph format for the UCSC genome browser.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0.8'
          id: break_point_threshold
          type: float?
          label: Break point threshold
          doc: >-
            Positive value of threshold for segmentation of normalized profiles.
            The closer it is to zero, the more breakpoints will be called. Its
            recommended value is between 0.1 and 1.2.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '2'
          id: break_point_type
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - '1'
                - '2'
                - '3'
                - '4'
              name: break_point_type
          label: Break point type
          doc: >-
            Desired behavior in the ambiguous regions (poly-N or low mappability
            regions between two different copy number values). 0: the "unknown"
            region is attached to the "known" region on the right  1: make a
            separate fragment of this “unknown” region and then attaches it to
            the left or to the right region choosing the longer one  2: make a
            separate fragment of this “unknown” region and then attaches it to
            the left or to the right region but the “ploidy” copy number has a
            priority  3: make a separate fragment of this “unknown” region and
            then attaches it to the left or to the right region choosing the
            longer one but this “known” region should make at least half-size of
            the “unknown” region  4: make a separate fragment of this “unknown”
            region and do not assign any copy number to this region at all
        - 'sbg:category': Target
          id: capture_regions
          type: File?
          label: Capture regions
          doc: >-
            Capture regions in .bed format; sorted .bed file should contain the
            following colomns: chr   0-based start   1-based end.
          'sbg:fileTypes': BED
        - 'sbg:category': File Input
          id: chr_len
          type: File
          label: Chromosomes length file
          doc: Chromosome length file in a tab-delimited format.
          'sbg:fileTypes': 'TXT, LEN, SIZES'
        - 'sbg:toolDefaultValue': '0.05'
          id: coeff_var
          type: float?
          label: Coefficient of variation
          doc: Coefficient of variation to evaluate necessary window size.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0'
          id: contamination
          type: float?
          label: Contamination
          doc: >-
            A priori known value of tumor sample contamiantion by normal cells.
            Set "contaminationAdjustment=TRUE" to correct for the contamination.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: contamination_adjustment
          type: boolean?
          label: Contamination adjustment
          doc: >-
            Set TRUE to correct for contamination by normal cells. If
            "contamination" is not provided, it will automatically evaluate the
            level of contamination.
        - 'sbg:category': General
          'sbg:toolDefaultValue': >-
            3&4 (GC-content based normalization, WGS) or 1
            (control-read-count-based normalization, WES)
          id: degree
          type: float?
          label: Degree
          doc: Degree of polynomial.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'WGS: 0, WES: 1'
          id: force_GC_content_normalization
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - '1'
                - '2'
              name: force_GC_content_normalization
          label: Force GC content normalization
          doc: >-
            Set 1 or 2 to correct the Read Count (RC) for GC-content bias and
            low mappability even when you have a control sample. 0: simply model
            "sample RC ~ Control RC"  1: normalize the sample and the control RC
            using GC-content and then calculate the ratio "Sample RC/contol RC"
            2: model "sample RC ~ Control RC" bias, and then normalize for
            GC-content.
        - 'sbg:category': File Input
          id: gem_mappability_file
          type: File?
          label: GEM mappability file
          doc: Mappability file in GEM format.
          'sbg:fileTypes': GEM
        - 'sbg:category': Execution
          'sbg:toolDefaultValue': '1 - with GC-content, 0 - with a control dataset'
          id: intercept
          type: float?
          label: Intercept
          doc: Intercept of polynomial.
        - 'sbg:category': Control
          id: mate_copynumber_file_control
          type: File?
          label: Mate copy number file control
          doc: >-
            Raw copy number profile for a given window-size (higher priority
            than mateFile)  (don't need to provide a mateFile if
            mateCopyNumberFile is provided).
          'sbg:fileTypes': CPN
        - 'sbg:category': Sample
          id: mate_copynumber_file_sample
          type: File?
          label: Mate copy number file sample
          doc: >-
            Raw copy number profile for a given window-size (higher priority
            than mateFile)  (don't need to provide a mateFile if
            mateCopyNumberFile is provided).
          'sbg:fileTypes': CPN
        - 'sbg:category': Control
          id: mate_file_control
          type: File?
          label: Mate file Control
          doc: >-
            Mapped reads (can be single end reads, mate-pairs or paired-end
            reads).
          'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
        - 'sbg:category': Sample
          id: mate_file_sample
          type: File?
          label: Mate file Sample
          doc: >-
            Mapped reads (can be single end reads, mate-pairs or paired-end
            reads).
          'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
        - 'sbg:category': Control
          id: mate_orientation_control
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - RF
                - FR
                - FF
              name: mate_orientation_control
          label: Mate orientation control
          doc: >-
            Format of reads (in mateFile). 0 (for single ends), RF (Illumina
            mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
        - 'sbg:category': Sample
          id: mate_orientation_sample
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - RF
                - FR
                - FF
              name: mate_orientation_sample
          label: Mate orientation sample
          doc: >-
            Format of reads (in mateFile). 0 (for single ends), RF (Illumina
            mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
        - 'sbg:category': General
          'sbg:toolDefaultValue': 0.55 (change only if you run Control-FREEC on a bacterial genome)
          id: max_expected_GC
          type: float?
          label: Maximum expected GC
          doc: >-
            Maximal exptected value of the GC-content for the prior evaluation
            of "Read Count ~ GC-content" dependancy.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '8'
          id: max_threads
          type: int?
          label: Maximum threads
          doc: Number of threads (multi-threading mode).
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'WES: 3, WGS: 1'
          id: min_CNA_length
          type: int?
          label: Minimum CNA length
          doc: Minimal number of consecutive windows to call a CNA.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 0.35 (change only if you run Control-FREEC on a bacterial genome)
          id: min_expected_GC
          type: float?
          label: Minimal exptected GC
          doc: >-
            Minimal expected value of the GC-content for the prior evaluation of
            "Read Count ~ GC-content" dependency.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0.85'
          id: min_map_per_w
          type: float?
          label: Minimum mappability per window
          doc: >-
            Only windows with fraction of mappable positions higher than or
            equal to this threshold will be considered  (if "gemMappabilityFile"
            is not provided, one uses the percentage of non-N letters per
            window).
        - 'sbg:category': General
          'sbg:toolDefaultValue': >-
            100 (meaning "do not look for subclones"). Suggested: 20 (or 0.2)
            for WGS and 30 (or 0.3) for WES.
          id: min_subclone_presence
          type: int?
          label: Minimal subclone presence
          doc: Detects subclones present in x% of cell population
        - 'sbg:category': Control
          id: mini_pileup_control
          type: File?
          label: Mini pileup Control
          doc: >-
            Mini pileup file created from the corresponding BAM file dring a
            previous run of Control-FREEC - providing this file will
            significantly speed up the whole process.
          'sbg:fileTypes': PILEUP
        - 'sbg:category': Sample
          id: mini_pileup_sample
          type: File?
          label: Mini pileup Sample
          doc: >-
            Mini pileup file created from the corresponding BAM file dring a
            previous run of Control-FREEC - providing this file will
            significantly speed up the whole process.
          'sbg:fileTypes': PILEUP
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: minimal_coverage_per_position
          type: int?
          label: Minimal coverage per position
          doc: >-
            Minimal read coverage for a position to be considered in BAF
            analysis.
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: minimal_quality_per_position
          type: int?
          label: Minimal quality per position
          doc: >-
            Minimal sequencing quality for a position to be considered in BAF
            analysis. Default: 0; using this option can slow down reading of
            pileup files.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: noisy_data
          type: boolean?
          label: Noisy data
          doc: >-
            Set TRUE for target resequencing data (e.g., exome-seq) to avoid
            false positive predictions due to nonuniform capture.
        - 'sbg:category': General
          id: ploidy
          type: 'int[]'
          label: Ploidy
          doc: >-
            Genome ploidy. In case of doubt, you can set different values and
            Control-FREEC will select the one that explains most observed CNAs
            (eg. 2,3,4)
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'TRUE'
          id: print_NA
          type: boolean?
          label: Print NA
          doc: >-
            Set FALSE to avoid printing "-1" to the _ratio.txt files Useful for
            exome-seq or targeted sequencing data.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 10; recommended value >=50 for for exome data
          id: read_cnt_threshold
          type: int?
          label: Read count threshold
          doc: >-
            Threshold on the minimal number of reads per window in the control
            sample Useful for exome-seq or targeted sequencing data.
        - 'sbg:category': BAF
          id: reference
          type: File
          label: Reference file
          doc: >-
            Reference file that will be divided as needed to separate
            chromosomes and contigs.
          'sbg:fileTypes': 'FA, FASTA'
        - 'sbg:category': General
          id: sex
          type:
            - 'null'
            - type: enum
              symbols:
                - XX
                - XY
              name: sex
          label: Sex
          doc: Sample sex.
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: shift_in_quality
          type: int?
          label: Shift in quality
          doc: >-
            Basis for Phred quality. Default: 0; usually 33 or 64; see fastq
            quality.
        - 'sbg:category': BAF
          id: snp_file
          type: File?
          label: Known SNPs
          doc: Known SNPs.
          'sbg:fileTypes': 'TXT, VCF'
        - 'sbg:category': General
          id: step
          type: int?
          label: Step
          doc: 'Step (used only when "window" is specified - Ex: 10000).'
        - 'sbg:category': General
          'sbg:toolDefaultValue': '50000'
          id: telocentromeric
          type: int?
          label: Telocentromeric
          doc: >-
            Length of pre-telomeric and pre-centromeric regions: Control-FREEC
            will not output small CNAs and LOH found within these regions (they
            are likely to be false because of mappability/genome assembly
            issues) 50000 is OK for human/mouse genomes. Use smaller values for
            yeasts and flies.
        - 'sbg:category': Execution
          id: total_memory
          type: int?
          label: 'Total memory [MB]'
          doc: Total amount of memory in MB reserved on the instance.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: unique_match
          type: boolean?
          label: Unique match
          doc: >-
            Use a mappability profile to correct read counts (in this case a
            mappability file must be provided with "gemMappabilityFile").
        - id: window
          type: int?
          label: Window
          doc: explicit window size (higher priority than coefficientOfVariation)
      outputs:
        - id: GC_profile
          doc: GC profile output file.
          label: GC profile
          type: File?
          outputBinding:
            glob: |-
              ${
                  return "GC_profile.targetedRegions.cnp"
              }
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': CNP
        - id: cnvs
          doc: File with coordinates of predicted copy number alterations.
          label: CNVs output
          type: File?
          outputBinding:
            glob: '*_CNVs'
            outputEval: |+
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }

          'sbg:fileTypes': TXT
        - id: cnvs_pvalue
          doc: >-
            File with coordinates of predicted copy number alterations with
            p-values.
          label: CNVs with p-value
          type: File?
          outputBinding:
            glob: '*.p.value.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: config_script
          doc: Configuration script used for running.
          label: Configuration script used for running
          type: File?
          outputBinding:
            glob: config.txt
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: control_cpn
          doc: >-
            Control CPN file to be used in the future for more efficient
            computation.
          label: Control CPN
          type: File?
          outputBinding:
            glob: '*_control.cpn'
            outputEval: |-
              ${  if (inputs.mate_file_control){
                  return inheritMetadata(self, inputs.mate_file_control)
              }
              }
          'sbg:fileTypes': CPN
        - id: control_pileup
          doc: >-
            Mini Pileup created for BAF calculation. It is used to speed up
            consequent runs with the same samples.
          label: Control Pileup
          type: File?
          outputBinding:
            glob: |-
              ${
                  if (inputs.mate_file_control) {

                      return inputs.mate_file_control.path.split('/').pop() + '_minipileup.pileup'

                  }
              }
            outputEval: |-
              ${
                  if (inputs.mate_file_control){
                  return inheritMetadata(self, inputs.mate_file_control)
              }
              }
          'sbg:fileTypes': PILEUP
        - id: info_txt
          doc: Parsable file with information about FREEC run.
          label: Info TXT
          type: File?
          outputBinding:
            glob: '*_info.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: pngs
          doc: >-
            Visalized normalized copy number profile with predicted CNAs as well
            as BAF profile (if dbSNP is provided)
          label: Copy number profile
          type: 'File[]?'
          outputBinding:
            glob: '*png'
          'sbg:fileTypes': PNG
        - id: ratio
          doc: >-
            File with ratios and predicted copy number alterations for each
            window.
          label: Ratio
          type: File?
          outputBinding:
            glob: '*_ratio.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: ratio_BedGraph
          doc: >-
            File with ratios in BedGraph format for visualization in the UCSC
            genome browser
          label: Ratio BedGraph
          type: File?
          outputBinding:
            glob: '*.BedGraph'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': BEDGRAPH
        - id: sample_BAF
          doc: >-
            File with B-allele frequencies for each possibly heterozygous SNP
            position
          label: BAF sample file
          type: File?
          outputBinding:
            glob: '*_BAF.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: sample_cpn
          doc: >-
            Sample CPN file to be used in the future for more efficient
            computation.
          label: Sample CPN
          type: File?
          outputBinding:
            glob: '*_sample.cpn'
            outputEval: |-
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': CPN
        - id: sample_pileup
          doc: >-
            Mini Pileup created for BAF calculation. It is used to speed up
            consequent runs with the same samples.
          label: Sample Pileup
          type: File?
          outputBinding:
            glob: |-
              ${
                  if (inputs.mate_file_sample) {

                      return inputs.mate_file_sample.path.split('/').pop() + '_minipileup.pileup'
                  }
              }
            outputEval: |-
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': PILEUP
      doc: >-
        Control-FREEC analyzes copy-number variants and allelic imbalances in
        exome and whole-genome DNA sequencing.


        This tool automatically computes, normalizes and segments copy number
        and beta allele frequency (BAF) profiles, then calls copy number
        alterations and LOH. [1]


        *A list of **all inputs and parameters** with corresponding descriptions
        can be found at the bottom of the page.*

        ### Common Use Cases


        * The **chrLenFile** input is required and can be found in Public
        Reference Files as **Homo\_sapiens\_assembly38.fasta.sizes**,
        **ucsc.hg19.fasta.sizes** and **human\_g1k\_v37\_decoy.fasta.sizes**.


        * The **ploidy** parameter is required. In case of doubt, different
        values can be set and Control-FREEC will select the one that explains
        the most observed CNVs.


        * Normal and control sample can be provided through two possible inputs:
             * **mateFile**, a file with mapped reads
             * **mateCopyNumberFile**, a raw copy number file created for both normal and control sample, provided through **mateFile** in a first run, and can be reused in the future runs for more efficient computation.



        * **A control (matched normal) sample is optional for whole genome
        sequencing data but mandatory for whole exome or targeted sequencing
        data.**


        * Similar to **mateCopyNumberFile**, a **Mini pileup Sample** and **Mini
        pileup Control** files can be created in the first run, if the **Known
        SNPs** file is provided. Consequently, by providing these files as
        inputs in future tasks, execution time will decrease significantly.


        * If a **mateFile** is specified, the **mateOrientation** parameter must
        be set.


        * In order to create a **BAF profile**, one of the following options
        must be implemented:
            * **mateFile** + **Known SNPs** 
            * **mateCopyNumberFile** + **mateFile** + **KnownSNPs**
            * **mateCopyNumberFile** + **miniPileup** + **KnownSNPs**

        ### Changes Introduced by Seven Bridges


        * Based on the input parameters, a config file is created in order to
        properly run Control-FREEC.


        ### Common Issues and Important Notes


        * **A control (matched normal) sample is optional for whole genome
        sequencing data but mandatory for whole exome or targeted sequencing
        data.**


        * A **gemMappabilityFile** can be used only in the mode without a
        control sample.


        * If a **mateFile** is specified, the **mateOrientation** parameter must
        be set.


        * Currently, there is an issue with creating a **BAF sample file** with
        the b37 notation. The genotypes for CNV regions are, however, created.



        ### Performance Benchmarking


        The instance set for this tool is the AWS c4.2xlarge instance with 8
        vCPUs, 15 GiB of RAM and 1 TB of EBS (disk space).

        |     BAM size in GB    | Type |  Instance  | Duration | Cost ($) |

        |:--------------------:|:----:|:----------:|:--------:|:--------:|

        |         2x12 (Normal-Tumor)         |  WES | c4.2xlarge |  1h 52m 
        |    0.8   |

        |   100 (Tumor-only)   |  WGS | c4.2xlarge |  17h 43m |     7    |

        | 2x100 (Normal-Tumor) |  WGS | c4.2xlarge |   1d 8h  |    13    |

        |   100 (Tumor-only)   |  WGS | c4.8xlarge |  6h 30m |     10    |

        | 2x100 (Normal-Tumor) |  WGS | c4.8xlarge |   11h  |    18    |


        An instance with more resources can be obtained by providing inputs for
        **Maximum threads** and **Total memory [MB]**.


        *Cost can be significantly reduced by using **spot instances**. Visit
        the [Knowledge
        Center](https://docs.sevenbridges.com/docs/about-spot-instances) for
        more details.*  


        ###References

        [1] [Control-FREEC: Prediction of copy number alterations and loss of
        heterozygosity using deep-sequencing
        data](http://boevalab.com/FREEC/tutorial.html#install)
      label: Control-FREEC 11.6
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            ${
                // script for splitting the genome fasta into chromosomes
                return 'python split_fasta.py ' + inputs.reference.path

            }
        - position: 1
          shellQuote: false
          valueFrom: '&&'
        - position: 2
          prefix: ''
          shellQuote: false
          valueFrom: /opt/controlfreec/FREEC/src/freec
        - position: 3
          shellQuote: false
          valueFrom: '-conf'
        - position: 4
          shellQuote: false
          valueFrom: config.txt
        - position: 5
          shellQuote: false
          valueFrom: '&&'
        - position: 6
          shellQuote: false
          valueFrom: |-
            ${




                if (inputs.mate_file_sample) {
                    filepath = inputs.mate_file_sample.path
                    filename = filepath.split("/").pop()
                } else {
                    filepath = inputs.mate_copynumber_file_sample.path
                    filename = filepath.split("/").pop()
                }

                CNVs = filename + "_CNVs"
                ratio = filename + "_ratio" + ".txt"


                return "cat assess_significance.R | R --slave --args " + CNVs + " " + ratio
            }
        - position: 7
          shellQuote: false
          valueFrom: '&&'
        - position: 8
          shellQuote: false
          valueFrom: |-
            ${
                return "line=$(cat *info.txt | grep Output_Ploidy | sed -E 's/.+([0-9]+)/\\1/')"
            }
        - position: 9
          shellQuote: false
          valueFrom: '&&'
        - position: 10
          shellQuote: false
          valueFrom: |-
            ${
                return "cat makeGraph.R | R --slave --args"
            }
        - position: 11
          shellQuote: false
          valueFrom: $line
        - position: 12
          shellQuote: false
          valueFrom: |-
            ${
                if (inputs.mate_file_sample) {
                    filepath = inputs.mate_file_sample.path
                    filename = filepath.split("/").pop()
                } else {
                    filepath = inputs.mate_copynumber_file_sample.path
                    filename = filepath.split("/").pop()
                }

                ratio = filename + "_ratio" + ".txt"

                return ratio
            }
        - position: 13
          shellQuote: false
          valueFrom: |-
            ${

                if (inputs.snp_file) {

                    sufix = "_BAF"
                    sufix_ext = ".txt"

                    if (inputs.mate_file_sample) {
                        filepath = inputs.mate_file_sample.path
                        filename = filepath.split("/").pop()
                    } else {
                        filepath = inputs.mate_copynumber_file_sample.path
                        filename = filepath.split("/").pop()
                    }


                    new_filename = filename + sufix + sufix_ext

                    return new_filename
                }
            }
        - position: 114
          shellQuote: false
          valueFrom: |-
            ${ //conversion of file names

                if (inputs.mate_file_control) {
                    if (inputs.mate_file_control.path.split('.').pop() != 'pileup') {
                        com = ''
                        com += '&& mv sample.pileup '

                    }
                }


            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: |-
            ${
                if (inputs.total_memory) {
                    return inputs.total_memory
                } else {
                    return 15000
                }
            }
          coresMin: |-
            ${
                if (inputs.max_threads) {
                    return inputs.max_threads
                } else {
                    return 8
                }
            }
        - class: DockerRequirement
          dockerImageId: caf6947244fa
          dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
        - class: InitialWorkDirRequirement
          listing:
            - entryname: config.txt
              entry: |-
                ${

                    // The function returns the concatenated line for config file
                    function makeline(content, p1, p2) {
                        if (p2 != null) {
                            if (p2.path != null) {
                                p2 = p2.path
                            }
                            content = content.concat(p1)
                            content = content.concat(" = ")
                            content = content.concat(p2)
                            content = content.concat("\n")
                        }
                        return content

                    }

                    // General section
                    content = "[general]\n\n"
                    content = makeline(content, "BedGraphOutput", inputs.bed_graph_output)
                    content = content.concat("bedtools = /opt/bedtools2/bin/bedtools\n")
                    content = makeline(content, "chrLenFile", inputs.chr_len)
                    content = makeline(content, "breakPointThreshold", inputs.break_point_threshold)
                    content = makeline(content, "breakPointType", inputs.break_point_type)
                    content = makeline(content, "chrFiles", ".")
                    content = makeline(content, "coefficientOfVariation", inputs.coeff_var)

                    if (inputs.capture_regions) {
                        content = content.concat("window = 0\n")
                    } else {
                        content = makeline(content, "window", inputs.window)
                    }

                    content = makeline(content, "contamination", inputs.contamination)
                    content = makeline(content, "contaminationAdjustment", inputs.contamination_adjustment)
                    content = makeline(content, "degree", inputs.degree)
                    content = makeline(content, "forceGCcontentNormalization", inputs.force_GC_content_normalization)
                    content = makeline(content, "GCcontentProfile", inputs.GC_content_profile)
                    content = makeline(content, "gemMappabilityFile", inputs.gem_mappability_file)
                    content = makeline(content, "intercept", inputs.intercept)
                    content = makeline(content, "minCNAlength", inputs.min_CNA_length)
                    content = makeline(content, "minMappabilityPerWindow", inputs.min_map_per_w)
                    content = makeline(content, "minExpectedGC", inputs.min_expected_GC)
                    content = makeline(content, "maxExpectedGC", inputs.max_expected_GC)
                    content = makeline(content, "minimalSubclonePresence", inputs.min_subclone_presence)
                    if (inputs.max_threads) {
                        content = makeline(content, "maxThreads", inputs.max_threads)
                        content = makeline(content, "SambambaThreads", inputs.max_threads)
                    } else {
                        content = content.concat("maxThreads = 8\n")
                        content = content.concat("SambambaThreads = 8\n")
                    }
                    content = makeline(content, "noisyData", inputs.noisy_data)
                    content = makeline(content, "ploidy", inputs.ploidy.toString())
                    content = makeline(content, "printNA", inputs.print_NA)
                    content = makeline(content, "readCountThreshold", inputs.read_cnt_threshold)
                    content = content.concat("sambamba = /opt/sambamba_0.5.9/sambamba_v0.5.9\n")

                    content = content.concat("samtools = /opt/samtools-1.3.1/samtools\n")
                    content = makeline(content, "sex", inputs.sex)
                    content = makeline(content, "step", inputs.step)
                    content = makeline(content, "telocentromeric", inputs.telocentromeric)
                    content = makeline(content, "uniqueMatch", inputs.unique_match)


                    // Sample section

                    content = content.concat("\n[sample]\n\n")
                    content = makeline(content, "mateFile", inputs.mate_file_sample)
                    content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_sample)
                    content = makeline(content, "miniPileup", inputs.mini_pileup_sample)
                    if (inputs.mate_file_sample) {
                        if (inputs.mate_file_sample.path.split('.').pop() == "gz") {
                            content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').slice(-2, -1)[0])
                        } else {
                            content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').pop())
                        }
                        content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
                    }


                    // Control section

                    content = content.concat("\n[control]\n\n")
                    content = makeline(content, "mateFile", inputs.mate_file_control)
                    content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_control)
                    content = makeline(content, "miniPileup", inputs.mini_pileup_control)
                    if (inputs.mate_file_control) {
                        content = makeline(content, "inputFormat", inputs.mate_file_control.path.split('.').pop())
                        content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
                    }




                    // BAF section

                    content = content.concat("\n[BAF]\n\n")
                    content = makeline(content, "minimalCoveragePerPosition", inputs.minimal_coverage_per_position)
                    content = makeline(content, "minimalQualityPerPosition", inputs.minimal_quality_per_position)
                    content = makeline(content, "shiftInQuality", inputs.shift_in_quality)
                    if (inputs.snp_file) {
                        content = makeline(content, "SNPfile", inputs.snp_file)
                        if (inputs.mate_file_sample) {
                            if ((inputs.mate_file_sample.path.split('.').pop().toUpperCase() != 'PILEUP') &&
                                (inputs.mate_file_sample.path.split('.').slice(-2, -1)[0].toUpperCase() != 'PILEUP')) {
                                content = makeline(content, "makePileup", inputs.snp_file)
                                content = makeline(content, "fastaFile", inputs.reference)
                            }
                        }
                    }

                    // Target section

                    content = content.concat("\n[target]\n\n")
                    content = makeline(content, "captureRegions", inputs.capture_regions)

                    return content
                }
              writable: false
            - entryname: split_fasta.py
              entry: |-
                import sys

                with open(sys.argv[1], "r") as f:
                    fasta = f.readlines()

                ref_lines = {}
                for i in range(0, len(fasta)):
                    if fasta[i][0] == ">":
                        chrom = fasta[i].split()[0].split(">")[1]
                        print("Reading chromosome: " + chrom)
                        ref_lines[chrom] = [fasta[i]]
                    else:
                        ref_lines[chrom].append(fasta[i])

                for chromosome, lines in ref_lines.items():
                    print("Creating " + chromosome + ".fasta")
                    with open(chromosome + ".fasta", "w") as chr_fasta:
                        for line in lines:
                            chr_fasta.write(line)
              writable: false
            - entryname: assess_significance.R
              entry: "#!/usr/bin/env Rscript\n\nlibrary(rtracklayer)\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\nratio<-data.frame(dataTable)\n\ndataTable <- read.table(args[4], header=FALSE)\ncnvs<- data.frame(dataTable) \n\nratio$Ratio[which(ratio$Ratio==-1)]=NA\n\ncnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  \nratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)\n\noverlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)\nnormals <- setdiff(ratio.bed,cnvs.bed)\nnormals <- subsetByOverlaps(ratio.bed,normals)\n\n#mu <- mean(score(normals),na.rm=TRUE)\n#sigma<- sd(score(normals),na.rm=TRUE)\n\n#hist(score(normals),n=500,xlim=c(0,2))\n#hist(log(score(normals)),n=500,xlim=c(-1,1))\n\n#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))\n#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)\n\n#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))\n#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)\n\nnumberOfCol=length(cnvs)\n\nfor (i in c(1:length(cnvs[,1]))) {\n  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))\n  #wilcox.test(values,mu=mu)\n  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)\n\tif(class(resultw)==\"try-error\") return(list(\"statistic\"=NA,\"parameter\"=NA,\"p.value\"=NA,\"null.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultw}\n  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)\n\tif(class(resultks)==\"try-error\") return(list(\"statistic\"=NA,\"p.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultks}\n  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)\n  #\tif(class(resultks)==\"try-error\") NA) else resultks\n  cnvs[i,numberOfCol+1]=W(values,normals)$p.value\n  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value\n  }\n\nif (numberOfCol==5) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==7) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==9) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"somatic/germline\",\"precentageOfGermline\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nwrite.table(cnvs, file=paste(args[4],\".p.value.txt\",sep=\"\"),sep=\"\\t\",quote=F,row.names=F)"
              writable: false
            - entryname: makeGraph.R
              entry: "#!/usr/bin/env Rscript\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\n\nratio<-data.frame(dataTable)\nploidy <- type.convert(args[4])\n\n\npng(filename = paste(args[5],\".log2.png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile (log2)\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[136])\n\t\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy), pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\t#points(ratio$Start[tt],log2(ratio$MedianRatio[tt]), pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\npng(filename = paste(args[5],\".png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nmaxLevelToPlot <- 3\nfor (i in c(1:length(ratio$Ratio))) {\n\tif (ratio$Ratio[i]>maxLevelToPlot) {\n\t\tratio$Ratio[i]=maxLevelToPlot;\n\t}\n}\n\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136])\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)\t\n\tpoints(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136],cex=4)\n\t \n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\t#points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\n\n\nif (length(args)>5) {\n\tdataTable <-read.table(args[6], header=TRUE);\n\tBAF<-data.frame(dataTable)\n\n\tpng(filename = paste(args[6],\".png\",sep = \"\"), width = 1180, height = 1180,\n\t    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\n\tplot(1:10)\n\top <- par(mfrow = c(5,5))\n\n\tfor (i in c(1:22,'X','Y')) {\n\t    tt <- which(BAF$Chromosome==i)\n\t    if (length(tt)>0){\n\t\tlBAF <-BAF[tt,]\n\t\tplot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste (\"position, chr\",i),ylab = \"BAF\",pch = \".\",col = colors()[1])\n\n\t\ttt <- which(lBAF$A==0.5)\t\t\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[92])\n\t\ttt <- which(lBAF$A!=0.5 & lBAF$A>=0)\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[62])\n\t\ttt <- 1\n\t\tpres <- 1\n\n\t\tif (length(lBAF$A)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$A)-pres-1))) {\n\t\t\t\tif (lBAF$A[j]==lBAF$A[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$A[tt],pch = \".\",col = colors()[24],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$B[tt],pch = \".\",col = colors()[24],cex=4)\t\n\t\t}\n\n\t\ttt <- 1\n\t\tpres <- 1\n\t\tif (length(lBAF$FittedA)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$FittedA)-pres-1))) {\n\t\t\t\tif (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedA[tt],pch = \".\",col = colors()[463],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedB[tt],pch = \".\",col = colors()[463],cex=4)\t\n\t\t}\n\n\t   }\n\n\t}\n\tdev.off()\n\n}"
              writable: false
        - class: InlineJavascriptRequirement
          expressionLib:
            - |-
              var updateMetadata = function(file, key, value) {
                  file['metadata'][key] = value;
                  return file;
              };


              var setMetadata = function(file, metadata) {
                  if (!('metadata' in file)) {
                      file['metadata'] = {}
                  }
                  for (var key in metadata) {
                      file['metadata'][key] = metadata[key];
                  }
                  return file
              };

              var inheritMetadata = function(o1, o2) {
                  var commonMetadata = {};
                  if (!Array.isArray(o2)) {
                      o2 = [o2]
                  }
                  for (var i = 0; i < o2.length; i++) {
                      var example = o2[i]['metadata'];
                      for (var key in example) {
                          if (i == 0)
                              commonMetadata[key] = example[key];
                          else {
                              if (!(commonMetadata[key] == example[key])) {
                                  delete commonMetadata[key]
                              }
                          }
                      }
                  }
                  if (!Array.isArray(o1)) {
                      o1 = setMetadata(o1, commonMetadata)
                  } else {
                      for (var i = 0; i < o1.length; i++) {
                          o1[i] = setMetadata(o1[i], commonMetadata)
                      }
                  }
                  return o1;
              };

              var toArray = function(file) {
                  return [].concat(file);
              };

              var groupBy = function(files, key) {
                  var groupedFiles = [];
                  var tempDict = {};
                  for (var i = 0; i < files.length; i++) {
                      var value = files[i]['metadata'][key];
                      if (value in tempDict)
                          tempDict[value].push(files[i]);
                      else tempDict[value] = [files[i]];
                  }
                  for (var key in tempDict) {
                      groupedFiles.push(tempDict[key]);
                  }
                  return groupedFiles;
              };

              var orderBy = function(files, key, order) {
                  var compareFunction = function(a, b) {
                      if (a['metadata'][key].constructor === Number) {
                          return a['metadata'][key] - b['metadata'][key];
                      } else {
                          var nameA = a['metadata'][key].toUpperCase();
                          var nameB = b['metadata'][key].toUpperCase();
                          if (nameA < nameB) {
                              return -1;
                          }
                          if (nameA > nameB) {
                              return 1;
                          }
                          return 0;
                      }
                  };

                  files = files.sort(compareFunction);
                  if (order == undefined || order == "asc")
                      return files;
                  else
                      return files.reverse();
              };
      successCodes:
        - 0
      temporaryFailCodes:
        - 1
      'sbg:categories':
        - Copy-Number-Analysis
      'sbg:cmdPreview': >-
        python split_fasta.py pat/sbgpro/path/to/reference.ext &&
        /opt/controlfreec/FREEC-11.5/src/freec -conf config.txt && cat
        assess_significance.R | R --slave --args chr_19.noDup0.pileup.gz_CNVs
        chr_19.noDup0.pileup.gz_ratio.txt && line=$(cat *info.txt | grep
        Output_Ploidy | sed -E 's/.+([0-9]+)/\1/') && cat makeGraph.R | R
        --slave --args $line chr_19.noDup0.pileup.gz_ratio.txt
        chr_19.noDup0.pileup.gz_BAF.txt
      'sbg:image_url': null
      'sbg:license': GNU General Public License v3.0 only
      'sbg:links':
        - id: 'http://bioinfo-out.curie.fr/projects/freec/'
          label: Homepage
        - id: 'https://github.com/BoevaLab/FREEC'
          label: GitHub
        - id: >-
            http://bioinformatics.oxfordjournals.org/content/early/2011/12/05/bioinformatics.btr670
          label: Paper
      'sbg:toolAuthor': Bioinformatics Laboratory of Institut Curie
      'sbg:toolkit': Control-FREEC
      'sbg:toolkitVersion': '11.6'
      'sbg:projectName': BMS CNV v2 - Dev
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1560511086
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': ana_popic
          'sbg:modifiedOn': 1560511123
          'sbg:revisionNotes': ''
        - 'sbg:revision': 2
          'sbg:modifiedBy': vojislav_varjacic
          'sbg:modifiedOn': 1565953496
          'sbg:revisionNotes': updated docker image
        - 'sbg:revision': 3
          'sbg:modifiedBy': vojislav_varjacic
          'sbg:modifiedOn': 1566208044
          'sbg:revisionNotes': updated docker image to 11.6
        - 'sbg:revision': 4
          'sbg:modifiedBy': vojislav_varjacic
          'sbg:modifiedOn': 1566239550
          'sbg:revisionNotes': base command changed
      'sbg:appVersion':
        - v1.0
      'sbg:id': vtomic/bms-cnv-v2-dev/control-freec-11-5/4
      'sbg:revision': 4
      'sbg:revisionNotes': base command changed
      'sbg:modifiedOn': 1566239550
      'sbg:modifiedBy': vojislav_varjacic
      'sbg:createdOn': 1560511086
      'sbg:createdBy': ana_popic
      'sbg:project': vtomic/bms-cnv-v2-dev
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
      'sbg:contributors':
        - vojislav_varjacic
        - ana_popic
      'sbg:latestRevision': 4
      'sbg:publisher': sbg
      'sbg:content_hash': a16aedc2a87c562123a6bcc93b8ab88adbfcbfb2b33c713185baf17909209fac7
    label: Control-FREEC 11.6
    'sbg:x': -33
    'sbg:y': -222
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:projectName': MB ControlFREEC Troubleshoot
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1567606643
    'sbg:revisionNotes': Copy of bogdang/controlfreec/control-freec-loh-11-6/1
  - 'sbg:revision': 1
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1567607314
    'sbg:revisionNotes': added a couple outputs
  - 'sbg:revision': 2
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1567611916
    'sbg:revisionNotes': added allowance of addiotional instances
  - 'sbg:revision': 3
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1567612158
    'sbg:revisionNotes': updated secondary files for bam to match our style
  - 'sbg:revision': 4
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1567629364
    'sbg:revisionNotes': added more outputs for testing
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/brownm28/mb-controlfreec-troubleshoot/control-freec-loh-11-6/4.png
'sbg:appVersion':
  - v1.0
'sbg:id': brownm28/mb-controlfreec-troubleshoot/control-freec-loh-11-6/4
'sbg:revision': 4
'sbg:revisionNotes': added more outputs for testing
'sbg:modifiedOn': 1567629364
'sbg:modifiedBy': brownm28
'sbg:createdOn': 1567606643
'sbg:createdBy': brownm28
'sbg:project': brownm28/mb-controlfreec-troubleshoot
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - brownm28
'sbg:latestRevision': 4
'sbg:publisher': sbg
'sbg:content_hash': a915131592ebb403c906f04f38ae75a8cfac5d7f029d4ad1216e216b495afdef9
