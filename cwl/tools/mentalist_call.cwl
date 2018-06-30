class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: mentalist_call
baseCommand:
  - mentalist
  - call
inputs:
  - id: db
    type: File
    inputBinding:
      position: 0
      prefix: '--db'
  - id: sample_fastq_1
    type: File
    inputBinding:
      position: 0
      prefix: '-1'
  - id: sample_fastq_2
    type: File?
    inputBinding:
      position: 0
      prefix: '-2'
  - id: output_votes
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--output_votes'
  - id: output_special
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--output_special'
outputs:
  - id: mentalist_calls
    type: File?
    outputBinding:
      glob: sequence_types.tsv
label: mentalist_call
arguments:
  - position: 0
    prefix: '-o'
    valueFrom: sequence_types.tsv
requirements:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/mentalist:0.1.8--0'
