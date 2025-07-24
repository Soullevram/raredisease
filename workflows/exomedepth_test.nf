nextflow.enable.dsl=2

include { EXOMEDEPTH_SUBWORKFLOW } from '../subworkflows/local/exomedepth'

workflow {

  /*
   * Inputs defined via config or parameters:
   *  - aligned_bams: tuple(meta, path(bam))
   *  - bed_file: path to BED
   *  - reference_fasta: path to reference genome
   */

  EXOMEDEPTH_SUBWORKFLOW(
    aligned_bams:       aligned_bams,
    bed_file:           bed_file,
    reference_fasta:    reference_fasta
  )

  EXOMEDEPTH_SUBWORKFLOW.out.exomedepth_results.view()
}
