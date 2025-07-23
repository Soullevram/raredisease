process EXOMEDEPTH_CALL {

  tag "${meta.id}"

  input:
    tuple val(meta), path(bam)
    path bed_file

  output:
    tuple val(meta), path("${meta.id}.exomedepth.txt")

  script:
  """
  Rscript ExomeDepthScript.R \\
    --bam ${bam} \\
    --bed ${bed_file} \\
    --out ${meta.id}.exomedepth.txt
  """
}
