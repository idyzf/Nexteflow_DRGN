import static nextflow.Nextflow.*

import dto.AnalysisInputFactory
import dto.Logging
import logic.DragenParamsMapUpdater
import logic.PublicS3DataRetriever
import logic.SampleSpecificSettings

nextflow.enable.dsl = 2

// Default params
def non_dragen_params = [
    is_cloud: true,
    is_test: false,
    fastqs: [],
    bams: [],
    bais: [],
    crams: [],
    crais: [],
    oras: [],
    ref_tar: '',
    is_raw_downsampling: false,
    additional_args: '',
    additional_vc_args: '',
    additional_files: [],
    enable_dragen_reports: false,
    samples_per_node: 1,
    public_s3_base_url: null,
    enable_splitvc: null, // controlled also by legacy mode
    enable_splitvc_spot: null, // controlled also by legacy mode
    legacy_mode: true,
]

// Load the defaults into the params object if the key does not exist
params.putAll(non_dragen_params.findAll { !params.containsKey(it.key) })

include { DRAGEN_GERMLINE } from './modules/germline' addParams(
    enable_splitvc: params.enable_splitvc != null ? params.enable_splitvc : !params.getOrDefault("legacy_mode", true),
    enable_splitvc_spot:
        params.enable_splitvc_spot != null ? params.enable_splitvc_spot : !params.getOrDefault("legacy_mode", true)
)
include { INTERSECT_METHYLATION_REPORTS } from './modules/intersect-methylation-reports'


def data_retriever = new PublicS3DataRetriever(
    System.getenv('CLOUD_REGION'), params.public_s3_base_url
)

workflow GERMLINE_ICA {
    take:
        params_ch // type: Map<String, Object>
        sample_specific_settings // SampleSpecificSettings

    main:
        params_ch.view { Logging.log_object("GERMLINE_ICA input", it) }

        germline_analysis_input_ch = params_ch
            .combine(sample_specific_settings)
            .map { Map<String, Object> params_map, SampleSpecificSettings sample_settings ->
                // Preprocess
                params_map = DragenParamsMapUpdater.update_germline_wgs_params_map(params_map, data_retriever)
                // Create germline
                return AnalysisInputFactory.create_germline_analysis_input(
                    params_map,
                    non_dragen_params.collect { it.key },
                    data_retriever,
                    true,
                    sample_settings
                )
            }

        DRAGEN_GERMLINE(
            germline_analysis_input_ch,
            params_ch
                .filter { Map<String, Object> params_map -> params_map.enable_dragen_reports }
                .map { Map<String, Object> params_map ->
                    // 5base/germline_wgs is the folder structure for DRAGEN reports
                    // Needs --methylation-conversion illumina option to be set
                    params_map.enable_5_base_methylation ? '5base/germline_wgs' : 'germline_wgs'
                }
        )

        // Intersect methylation report when cytosine bed file or vc target bed file is provided
        // and methylation_generate_cytosine_report true
        // todo: need to replace enable_5_base_methylation to methylation_generate_cytosine_report
        cytosine_bed_ch = params_ch
            .filter { it.enable_5_base_methylation && (it.vc_target_bed || it.cytosine_bed) }
            .map { file(it.cytosine_bed ?: it.vc_target_bed) }
        INTERSECT_METHYLATION_REPORTS(DRAGEN_GERMLINE.out.dragen.collect(), cytosine_bed_ch, false)

    emit:
        dragen = DRAGEN_GERMLINE.out.dragen
        dragen_report = DRAGEN_GERMLINE.out.dragen_report
}

workflow {
    GERMLINE_ICA(Channel.of(params), Channel.of(SampleSpecificSettings.create_null_object()))