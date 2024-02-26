##########################################ADDED SOME COMMENTS and FEEDBACK in code file#####################################################
#' summix_network
#'
#' @description
#' Helper function to plot the network diagram of estimated substructure proportions and similarity between reference groups
#'
#' @param data A dataframe of the observed and reference allele frequencies for N genetic variants. See data formatting document at \href{https://github.com/hendriau/Summix}{https://github.com/hendriau/Summix} for more information.
#' @param sum_res The resulting data frame from the summix function
#' @param reference A character vector of the column names for the reference groups.
#' @param N_reference numeric vector of the sample sizes for each of the K reference groups.
#' @param reference_colors A character vector of length K that specifies the color each reference group node in the network plot. If not specified, this defaults to K random colors.
#' 
#' @return network diagram with nodes as estimated substructure proportions and edges as degree of similarity between the given node pair
#' 
#' @importFrom magrittr "%>%"
#' @importFrom BEDASSLE "calculate.all.pairwise.Fst"
#' @importFrom visNetwork "visNetwork" "visNodes" "visEdges" "visLayout" "visExport"
#' @importFrom scales "percent"
#' @importFrom randomcoloR "distinctColorPalette"
#' 
#' @export

summix_network <- function(data = data,
                           sum_res = sum_res,
                           reference = reference,
                           N_reference = N_reference,
                           reference_colors=reference_colors){
  
  
  detected_refs <- names(sum_res[5:length(sum_res)])
  N_reference = N_reference[which(reference %in% detected_refs)]
  reference = reference[which(reference %in% detected_refs)]
  
  
  if (identical(which(sum_res[1,5:length(sum_res)]<.005), integer(0))){
    non_zero_refs <- reference
    non_zero_N_refs <- N_reference
  } else{
    non_zero_refs <- reference[-c(which(sum_res[1,5:length(sum_res)]<.005))]
    non_zero_N_refs <- N_reference[-c(which(sum_res[1,5:length(sum_res)]<.005))]
  }
  
  
  alleleCounts <- t(as.matrix(sweep(data[1:dim(data)[1],non_zero_refs], MARGIN=2, 2*non_zero_N_refs, `*`)))
  alleleNumber <- t(as.matrix(sweep(matrix(1, dim(data)[1], length(non_zero_refs)), MARGIN = 2, 2*non_zero_N_refs, `*`)))
  rownames(alleleNumber) <- non_zero_refs
  
  
  fst_pure <- calculate.all.pairwise.Fst(allele.counts = alleleCounts,
                                         sample.sizes = alleleNumber)
  
  colnames(fst_pure) <- non_zero_refs
  rownames(fst_pure) <- non_zero_refs
  
  edges <- data.frame(row=rownames(fst_pure)[row(fst_pure)[upper.tri(fst_pure)]],
                      col=colnames(fst_pure)[col(fst_pure)[upper.tri(fst_pure)]],
                      corr=fst_pure[upper.tri(fst_pure)])
  colnames(edges) <- c("from", "to", "weight")
  
  if (length(reference_colors) != 1){
    reference_colors = reference_colors[which(reference %in% non_zero_refs)]
  }else{
    reference_colors = distinctColorPalette(length(non_zero_refs))
    print(reference_colors)
  }
  
  nodes <- data.frame(id = non_zero_refs,
                      proportions = as.numeric(sum_res[1, non_zero_refs]),
                      references = non_zero_refs,
                      color = reference_colors)
  
  nodes$label <- paste0(nodes$references, " ", percent(round(nodes$proportions,2)))
  
  nodes$size <- 90*nodes$proportions
  edges$width = -1*as.numeric(5*scale(edges$weight))
  
  sn <- visNetwork(nodes, edges, width = "100%") %>%
    visNodes(
      shape = "dot",
      shadow = T,
      size = 100,
      font = list(size = 14),
      color = list(
        border = "#013848",
        highlight = "#FF8000"),
    ) %>%
    visEdges(
      shadow = FALSE,
      smooth = FALSE,
      color = list(color = "lightgray", highlight = "#C62F4B")
    ) 
  return(visExport(sn, name = "Summix_Network", type = "png"))
  
}



