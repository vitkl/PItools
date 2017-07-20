##' Extract binding region feature from cleaned MITAB2.7
##' @name MITABregionFeature
##' @import data.table
##' @param mitab data.table containing data in MITAB2.7 format, standard column names should be made data.table-compatible
##' @author Vitalii Kleshchevnikov
##' @return data.table to which the following columns were added: binding_region_A_type, binding_region_A, binding_region_B_type, binding_region_B
MITABregionFeature = function(mitab){
  # extracting Features_interactor_A
  mitab = copy(mitab)
  mitab = mitab[, .(pair_id,
                    IDs_interactor_A, IDs_interactor_B,
                    interactor_IDs_databases_A, interactor_IDs_databases_B,
                    Taxid_interactor_A, Taxid_interactor_B,
                    Publication_Identifiers, Confidence_values,
                    Host_organisms,
                    bait_prey_status_A, bait_prey_status_B,
                    Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                    Features_interactor_A, Features_interactor_B,
                    Identification_method_participant_A, Identification_method_participant_B,
                    Features_interactor_A_temp = unlist(strsplit(Features_interactor_A, "\\|"))), by = 1:NROW(mitab)]
  # extracting Features_interactor_B
  mitab = mitab[, .(pair_id,
                    IDs_interactor_A, IDs_interactor_B,
                    interactor_IDs_databases_A, interactor_IDs_databases_B,
                    Taxid_interactor_A, Taxid_interactor_B,
                    Publication_Identifiers, Confidence_values,
                    Host_organisms,
                    bait_prey_status_A, bait_prey_status_B,
                    Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                    Features_interactor_A, Features_interactor_B,
                    Identification_method_participant_A, Identification_method_participant_B, Features_interactor_A_temp,
                    Features_interactor_B_temp = unlist(strsplit(Features_interactor_B, "\\|"))), by = 1:NROW(mitab)]
  # extracting binding regions from Features_interactor_A
  mitab[, binding_associated_region_interactor_A := gsub("^binding-associated region:|\\(.*\\)$", "", grep("binding-associated region:", Features_interactor_A_temp, value = T)), by = Features_interactor_A_temp]
  mitab[, necessary_binding_region_interactor_A := gsub("^necessary binding region:|\\(.*\\)$", "", grep("necessary binding region:", Features_interactor_A_temp, value = T)), by = Features_interactor_A_temp]
  mitab[, sufficient_binding_region_interactor_A := gsub("^sufficient binding region:|\\(.*\\)$", "", grep("sufficient binding region:", Features_interactor_A_temp, value = T)), by = Features_interactor_A_temp]
  mitab[, direct_binding_region_interactor_A := gsub("^direct binding region:|\\(.*\\)$", "", grep("direct binding region:", Features_interactor_A_temp, value = T)), by = Features_interactor_A_temp]
  # collapsing binding regions from Features_interactor_A into 1 column keeping the most detailed annotation
  mitab[!is.na(direct_binding_region_interactor_A), binding_region_A_type := "direct_binding_region"]
  mitab[!is.na(direct_binding_region_interactor_A), binding_region_A := direct_binding_region_interactor_A]
  mitab[!is.na(sufficient_binding_region_interactor_A) & is.na(binding_region_A), binding_region_A_type := "sufficient_binding_region"]
  mitab[!is.na(sufficient_binding_region_interactor_A) & is.na(binding_region_A), binding_region_A := sufficient_binding_region_interactor_A]
  mitab[!is.na(necessary_binding_region_interactor_A) & is.na(binding_region_A), binding_region_A_type := "necessary_binding_region"]
  mitab[!is.na(necessary_binding_region_interactor_A) & is.na(binding_region_A), binding_region_A := necessary_binding_region_interactor_A]
  mitab[!is.na(binding_associated_region_interactor_A) & is.na(binding_region_A), binding_region_A_type := "binding_associated_region"]
  mitab[!is.na(binding_associated_region_interactor_A) & is.na(binding_region_A), binding_region_A := binding_associated_region_interactor_A]

  # extracting binding regions from Features_interactor_B
  mitab[, binding_associated_region_interactor_B := gsub("^binding-associated region:|\\(.*\\)$", "", grep("binding-associated region:", Features_interactor_B_temp, value = T)), by = Features_interactor_B_temp]
  mitab[, necessary_binding_region_interactor_B := gsub("^necessary binding region:|\\(.*\\)$", "", grep("necessary binding region:", Features_interactor_B_temp, value = T)), by = Features_interactor_B_temp]
  mitab[, sufficient_binding_region_interactor_B := gsub("^sufficient binding region:|\\(.*\\)$", "", grep("sufficient binding region:", Features_interactor_B_temp, value = T)), by = Features_interactor_B_temp]
  mitab[, direct_binding_region_interactor_B := gsub("^direct binding region:|\\(.*\\)$", "", grep("direct binding region:", Features_interactor_B_temp, value = T)), by = Features_interactor_B_temp]
  # collapsing binding regions from Features_interactor_B into 1 column keeping the most detailed annotation
  mitab[!is.na(direct_binding_region_interactor_B), binding_region_B_type := "direct_binding_region"]
  mitab[!is.na(direct_binding_region_interactor_B), binding_region_B := direct_binding_region_interactor_B]
  mitab[!is.na(sufficient_binding_region_interactor_B) & is.na(binding_region_B), binding_region_B_type := "sufficient_binding_region"]
  mitab[!is.na(sufficient_binding_region_interactor_B) & is.na(binding_region_B), binding_region_B := sufficient_binding_region_interactor_B]
  mitab[!is.na(necessary_binding_region_interactor_B) & is.na(binding_region_B), binding_region_B_type := "necessary_binding_region"]
  mitab[!is.na(necessary_binding_region_interactor_B) & is.na(binding_region_B), binding_region_B := necessary_binding_region_interactor_B]
  mitab[!is.na(binding_associated_region_interactor_B) & is.na(binding_region_B), binding_region_B_type := "binding_associated_region"]
  mitab[!is.na(binding_associated_region_interactor_B) & is.na(binding_region_B), binding_region_B := binding_associated_region_interactor_B]

  # split multiple sufficient_binding_region and necessary_binding_region into necessary_binding_regions
  # spliting necessary_binding_regions Features_interactor_A
  mitab = mitab[, .(pair_id,
                    IDs_interactor_A, IDs_interactor_B,
                    interactor_IDs_databases_A, interactor_IDs_databases_B,
                    Taxid_interactor_A, Taxid_interactor_B,
                    Publication_Identifiers, Confidence_values,
                    Host_organisms,
                    bait_prey_status_A, bait_prey_status_B,
                    Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                    Features_interactor_A, Features_interactor_B,
                    Identification_method_participant_A, Identification_method_participant_B,
                    Features_interactor_A_temp, Features_interactor_B_temp,
                    binding_region_A, binding_region_B,
                    binding_region_A_type, binding_region_B_type,
                    binding_region_A_temp = unlist(strsplit(binding_region_A, ","))), by = 1:NROW(mitab)]
  # spliting necessary_binding_regions Features_interactor_B
  mitab = mitab[, .(pair_id,
                    IDs_interactor_A, IDs_interactor_B,
                    interactor_IDs_databases_A, interactor_IDs_databases_B,
                    Taxid_interactor_A, Taxid_interactor_B,
                    Publication_Identifiers, Confidence_values,
                    Host_organisms,
                    bait_prey_status_A, bait_prey_status_B,
                    Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                    Features_interactor_A, Features_interactor_B,
                    Identification_method_participant_A, Identification_method_participant_B,
                    Features_interactor_A_temp, Features_interactor_B_temp,
                    binding_region_A, binding_region_B,
                    binding_region_A_type, binding_region_B_type,
                    binding_region_A_temp,
                    binding_region_B_temp = unlist(strsplit(binding_region_B, ","))), by = 1:NROW(mitab)]
  # record split region status
  mitab[grep(",", binding_region_A), binding_region_A_type := paste0("split_", binding_region_A_type)]
  mitab[grep(",", binding_region_A), binding_region_A := binding_region_A_temp]
  mitab[grep(",", binding_region_B), binding_region_B_type := paste0("split_", binding_region_B_type)]
  mitab[grep(",", binding_region_B), binding_region_B := binding_region_B_temp]

  # delete "?-?" or "?-423" or "453-?" regions
  mitab[grep("\\?", binding_region_A), binding_region_A_type := "undetermined_position"]
  mitab[grep("\\?", binding_region_A), binding_region_A := NA]
  mitab[grep("\\?", binding_region_B), binding_region_B_type := "undetermined_position"]
  mitab[grep("\\?", binding_region_B), binding_region_B := NA]
  # delete n-terminal region information
  mitab[grep("n", binding_region_A), binding_region_A_type := "n_terminal_range"]
  mitab[grep("n", binding_region_A), binding_region_A := NA]
  mitab[grep("n", binding_region_B), binding_region_B_type := "n_terminal_range"]
  mitab[grep("n", binding_region_B), binding_region_B := NA]
  # delete c-terminal region information
  mitab[grep("c", binding_region_A), binding_region_A_type := "c_terminal_range"]
  mitab[grep("c", binding_region_A), binding_region_A := NA]
  mitab[grep("c", binding_region_B), binding_region_B_type := "c_terminal_range"]
  mitab[grep("c", binding_region_B), binding_region_B := NA]
  # delete greater or smaller range information
  mitab[grep(">|<", binding_region_A), binding_region_A_type := "greater_or_smaller_range"]
  mitab[, binding_region_A := gsub(">|<", "",binding_region_A)]
  mitab[grep(">|<", binding_region_B), binding_region_B_type := "greater_or_smaller_range"]
  mitab[, binding_region_B := gsub(">|<", "",binding_region_B)]
  # modify blurry ends ("75..75-88..88") to be strict ("75-88")
  mitab[grep("\\.\\.[[:digit:]]+-[[:digit:]]+\\.\\.", binding_region_A), binding_region_A_type := paste0("blurry_end_",binding_region_A_type)]
  mitab[grep("\\.\\.[[:digit:]]+-[[:digit:]]+\\.\\.", binding_region_B), binding_region_B_type := paste0("blurry_end_",binding_region_B_type)]
  mitab[, binding_region_A := gsub("\\.\\.[[:digit:]]+-[[:digit:]]+\\.\\.", "-", binding_region_A)]
  mitab[, binding_region_B := gsub("\\.\\.[[:digit:]]+-[[:digit:]]+\\.\\.", "-", binding_region_B)]
  # modify blurry ends ("75-88..88") to be strict ("75-88")
  mitab[grep("-[[:digit:]]+\\.\\.", binding_region_A), binding_region_A_type := paste0("blurry_end_end_",binding_region_A_type)]
  mitab[grep("-[[:digit:]]+\\.\\.", binding_region_B), binding_region_B_type := paste0("blurry_end_end_",binding_region_B_type)]
  mitab[, binding_region_A := gsub("-[[:digit:]]+\\.\\.", "-", binding_region_A)]
  mitab[, binding_region_B := gsub("-[[:digit:]]+\\.\\.", "-", binding_region_B)]
  # modify blurry ends ("75..75-88") to be strict ("75-88")
  mitab[grep("\\.\\.[[:digit:]]+-", binding_region_A), binding_region_A_type := paste0("blurry_end_start_",binding_region_A_type)]
  mitab[grep("\\.\\.[[:digit:]]+-", binding_region_B), binding_region_B_type := paste0("blurry_end_start_",binding_region_B_type)]
  mitab[, binding_region_A := gsub("\\.\\.[[:digit:]]+-", "-", binding_region_A)]
  mitab[, binding_region_B := gsub("\\.\\.[[:digit:]]+-", "-", binding_region_B)]

  # split region into the start and the end positions
  mitab[, c("binding_region_A_start", "binding_region_A_end") := tstrsplit(binding_region_A, "-")]
  mitab[, c("binding_region_B_start", "binding_region_B_end") := tstrsplit(binding_region_B, "-")]

  # select final set of columns
  mitab = mitab[, .(pair_id,
                    IDs_interactor_A, IDs_interactor_B,
                    interactor_IDs_databases_A, interactor_IDs_databases_B,
                    Taxid_interactor_A, Taxid_interactor_B,
                    Publication_Identifiers, Confidence_values,
                    Host_organisms,
                    bait_prey_status_A, bait_prey_status_B,
                    Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                    Features_interactor_A, Features_interactor_B,
                    Identification_method_participant_A, Identification_method_participant_B,
                    binding_region_A_start, binding_region_A_end, binding_region_B_start, binding_region_B_end,
                    binding_region_A_type, binding_region_B_type)]
  return(unique(mitab))
}
