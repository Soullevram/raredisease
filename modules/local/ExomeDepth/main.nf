process EXOMEDEPTH_CALL {

  tag "${meta.id}"

  input:
    tuple val(meta), path(test_bam)
    val(control_bams) // list of paths (strings)
    path bed_file
    path fasta

  output:
    tuple val(meta), path("${meta.id}.exomedepth.txt")

  script:
  // Convert control list to comma-separated string
  def controls_str = control_bams.join(',')
  """
  Rscript ExomeDepthScript.R \\
    --bam ${test_bam} \\
    --controls ${controls_str} \\
    --bed ${bed_file} \\
    --fasta ${fasta} \\
    --out ${meta.id}.exomedepth.txt
  """
}
