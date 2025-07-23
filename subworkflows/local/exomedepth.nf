nextflow.enable.dsl=2

include { EXOMEDEPTH_CALL } from '../../modules/local/ExomeDepth'

workflow EXOMEDEPTH_SUBWORKFLOW {

  take:
    aligned_bams   // tuple(meta, bam)
    bed_file
    reference_fasta

  main:

    aligned_bams
      .collect()
      .flatMap { all_samples ->
        all_samples.collect { test_sample ->
          def test_meta = test_sample[0]
          def test_bam  = test_sample[1]
          def control_bams = all_samples
                              .findAll { it[0].id != test_meta.id }
                              .collect { it[1].toString() } // string paths

          return tuple(test_meta, test_bam, control_bams, bed_file, reference_fasta)
        }
      } | EXOMEDEPTH_CALL

  emit:
    exomedepth_results
}
