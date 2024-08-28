//
// Subworkflow with functionality specific to the nf-core/phyloplace pipeline
//

import groovy.json.JsonSlurper

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --genome GRCh37 --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Function to validate channels from input samplesheet
//
/** Leaving this for the moment, keeping my workflow code
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect{ it.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
**/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    if ( ! params.phyloplace_input && ! ( params.phylosearch_input && params.search_fasta ) && ! ( params.id && params.queryseqfile && params.refseqfile && params.refphylogeny && params.model ) ) {
        error("Please provide either --input sample.csv or separatate parameters for input, see documentation or --help")
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "EPA-NG (Barbera et al. 2019",
            "Gappa (Czech et al. 2020",
            "HMMER (Eddy 2011)",
            "MultiQC (Ewels et al. 2016)",
            "MAFFT (Katoh et al. 2002",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Barbera, Pierre, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, and Alexandros Stamatakis. “EPA-Ng: Massively Parallel Evolutionary Placement of Genetic Sequences.” Systematic Biology 68, no. 2 (March 1, 2019): 365–69. https://doi.org/10.1093/sysbio/syy054.</li>",
            "<li>Czech, Lucas, Pierre Barbera, and Alexandros Stamatakis. “Genesis and Gappa: Processing, Analyzing and Visualizing Phylogenetic (Placement) Data.” Bioinformatics 36, no. 10 (May 1, 2020): 3263–65. https://doi.org/10.1093/bioinformatics/btaa070.</li>",
            "<li>Eddy, Sean R. “Accelerated Profile HMM Searches.” PLoS Comput Biol 7, no. 10 (October 20, 2011): e1002195. https://doi.org/10.1371/journal.pcbi.1002195.</li>",
            "<li>Katoh, Kazutaka, Kazuharu Misawa, Kei‐ichi Kuma, and Takashi Miyata. “MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast Fourier Transform.” Nucleic Acids Research 30, no. 14 (July 15, 2002): 3059–66. https://doi.org/10.1093/nar/gkf436.</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()
    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}
