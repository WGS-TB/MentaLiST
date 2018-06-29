class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: mentalist_download_pubmlst
baseCommand:
  - mentalist
  - download_pubmlst
inputs:
  - 'sbg:toolDefaultValue': '31'
    id: k
    type: int
    inputBinding:
      position: 0
      prefix: '-k'
  - id: scheme
    type: string
    inputBinding:
      position: 0
      prefix: '-s'
outputs:
  - id: mentalist_kmer_db
    type: File
    outputBinding:
      glob: mentalist.db
  - id: scheme_fasta_files
    type: File?
    outputBinding:
      glob: scheme_fasta_files
  - id: profile
    type: File?
    outputBinding:
      glob: mentalist.db.profile
label: mentalist_download_pubmlst
arguments:
  - position: 0
    prefix: '--threads'
    valueFrom: $(runtime.cores)
  - position: 0
    prefix: '-o'
    valueFrom: scheme_fasta_files
  - position: 0
    prefix: '--db'
    valueFrom: mentalist.db
requirements:
  - class: ResourceRequirement
    coresMin: 0
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/mentalist:0.1.8--0'
  - class: InlineJavascriptRequirement
