# Script to create cptac-ccrcc-ptm-adjusted.csv and cptac-ccrcc-protein.csv
# in the extdata folder, for the PTM-Analysis vignette.
#
# Input: FragPipe output of the CPTAC clear cell renal cell carcinoma
# (CCRCC) cohort, TMT, with a phosphopeptide-enriched experiment (ptm/) and a
# global profiling experiment (proteome/), each with an msstats.csv and an
# MSstatsTMT_annotation.csv. The CPTAC pan-cancer data are described in
# Li Y, et al. Proteogenomic data and resources for pan-cancer analysis.
# Cancer Cell 2023;41(8):1397-1406. doi:10.1016/j.ccell.2023.06.009
# and can be downloaded from https://pdc.cancer.gov/pdc/cptac-pancancer.
# The MSstatsPTM steps are those of the May Institute PTM case study
# (Wed_AM-ptm_case_study.Rmd). Set data_dir to the folder holding ptm/ and
# proteome/.
#
# Output: the "T vs N" (tumor vs normal) comparison of groupComparisonPTM(),
# without contaminant (contam_) and decoy (rev_) proteins, with the columns
# the vignette uses. Protein IDs keep FragPipe's form,
# sp|<accession>|<mnemonic>, followed by the site for PTMs.
library(MSstatsPTM)
data_dir = "CPTAC"

ptm_input = read.csv(file.path(data_dir, "ptm/msstats.csv"))
ptm_annotation = read.csv(file.path(data_dir, "ptm/MSstatsTMT_annotation.csv"))
protein_input = read.csv(file.path(data_dir, "proteome/msstats.csv"))
protein_annotation = read.csv(file.path(data_dir,
                                        "proteome/MSstatsTMT_annotation.csv"))

msstatsptm_input = FragPipetoMSstatsPTMFormat(
    ptm_input,
    annotation = ptm_annotation,
    input_protein = protein_input,
    annotation_protein = protein_annotation,
    use_unmod_peptides = FALSE,
    label_type = "TMT",
    protein_id_col = "Protein",
    peptide_id_col = "Peptide.Sequence",
    mod_id_col = "STY",
    localization_cutoff = 0.75,
    remove_unlocalized_peptides = TRUE,
    Purity_cutoff = 0.6,
    PeptideProphet_prob_cutoff = 0.7,
    useUniquePeptide = TRUE,
    rmPSM_withfewMea_withinRun = FALSE,
    rmPeptide_OxidationM = TRUE,
    rmProtein_with1Feature = FALSE,
    summaryforMultipleRows = sum
)
summarized = dataSummarizationPTM_TMT(
    msstatsptm_input,
    method = "msstats",
    global_norm = TRUE,
    global_norm.PTM = TRUE,
    reference_norm = TRUE,
    reference_norm.PTM = TRUE,
    remove_norm_channel = TRUE,
    remove_empty_channel = TRUE,
    MBimpute = TRUE,
    MBimpute.PTM = TRUE,
    maxQuantileforCensored = NULL
)
model = groupComparisonPTM(
    summarized,
    data.type = "TMT",
    contrast.matrix = "pairwise",
    moderated = FALSE,
    adj.method = "BH",
    log_base = 2
)

keep_tumor_vs_normal = function(results) {
    results = as.data.frame(results)
    results$Protein = as.character(results$Protein)
    results[results$Label == "T vs N" &
                !grepl("^(contam_|rev_)", results$Protein), ]
}
adjusted = keep_tumor_vs_normal(model$ADJUSTED.Model)
protein = keep_tumor_vs_normal(model$PROTEIN.Model)

adjusted = data.frame(Protein = adjusted$Protein, Label = adjusted$Label,
                      log2FC = signif(adjusted$log2FC, 4),
                      adj.pvalue = signif(adjusted$adj.pvalue, 4),
                      GlobalProtein = as.character(adjusted$GlobalProtein),
                      Adjusted = adjusted$Adjusted)
protein = data.frame(Protein = protein$Protein, Label = protein$Label,
                     log2FC = signif(protein$log2FC, 4),
                     adj.pvalue = signif(protein$adj.pvalue, 4))

write.csv(adjusted, "inst/extdata/cptac-ccrcc-ptm-adjusted.csv",
          row.names = FALSE)
write.csv(protein, "inst/extdata/cptac-ccrcc-protein.csv", row.names = FALSE)
