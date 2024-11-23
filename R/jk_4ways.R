#' JK Fourway Dataframe
#'
#' @param dataframe1 first dataframe, will correspond to "X" axis of four-way
#' @param dataframe2 second dataframe, will correspond to "Y" axis of four-way
#' @param fc_col column with fold change info to be plotted, usually Log2-FC or Log-FC
#' @param fdr_col column with p-value/q-value/FDR value info to be plotted
#' @param fc_thres threshold of fold change to consider. Default is 0 because preference is to still plot genes below the FC threshold as gray.
#' @param fdr_thres threshold of p/q/FDR value to plot, default 0.05
#' @return A data frame with shared significant genes/proteins across two contrasts, ready for plotting 4-ways
#' @export
#' @importFrom ggplot2 aes geom_point geom_hline geom_vline element_text scale_y_continuous scale_x_continuous expansion geom_segment annotate xlab
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% count
#' @examples jk_sharedDE(testvol, fc_col = "log2FC", fdr_col = "qVal") #basic volcano plot using dataframe "testvol"
#' @examples jk_sharedDE(testvol, fc_col = "log2FC", fdr_col = "qVal", genenames = FALSE) #basic volcano plot with unlabeled points
#'
#'
jk_sharedDE <- function (dataframe1,
                         dataframe2,
                         gene_col = "Gene",
                         fc_col = "Log2FC", #fold change column
                         fdr_col = "FDR", #p or q value column
                         fdr_thres = 0.05,
                         fc_thres = 0) {


#code developed for Log2FC and FDR, so renaming everything as such; in principle works with Log10FC, or p-val instead of FDR
names(dataframe1)[names(dataframe1) == {{fc_col}}] <- "Log2FC"
names(dataframe1)[names(dataframe1) == {{fdr_col}}] <- "FDR"
names(dataframe1)[names(dataframe1) == {{gene_col}}] <- "Gene"

names(dataframe2)[names(dataframe2) == {{fc_col}}] <- "Log2FC"
names(dataframe2)[names(dataframe2) == {{fdr_col}}] <- "FDR"
names(dataframe2)[names(dataframe2) == {{gene_col}}] <- "Gene"


dataframe1 %>%
  group_by(Gene) %>%
  count() %>% filter(n>1) %>%
  dplyr::select(Gene) %>% unlist() %>%
  {.->>multiples.df1}

dataframe2 %>%
  group_by(Gene) %>%
  count() %>% filter(n>1) %>%
  dplyr::select(Gene) %>% unlist() %>%
  {.->>multiples.df2}

if (length(multiples.df1) > 0){
  print("Warning, dataframe 1 has duplicate genes. Duplicates will be filtered in four-way dataframe.")
}

if (length(multiples.df2) > 0){
  print("Warning, dataframe 2 has duplicate genes. Duplicates will be filtered in four-way dataframe.")
}

multiples <- unique(c(multiples.df1, multiples.df2))

dataframe1 %>%
  filter(!Gene %in% multiples) %>%
  filter(abs(Log2FC) > fc_thres,
         FDR < fdr_thres) %>%
  dplyr::select(Gene) %>% unlist() %>% {.->>sig_prot_a}


dataframe2 %>%
  filter(!Gene %in% multiples) %>%
  filter(abs(Log2FC) > fc_thres,
         FDR < fdr_thres) %>%
  dplyr::select(Gene) %>% unlist() %>% {.->>sig_prot_b}

unique(c(sig_prot_a, sig_prot_b)) %>% {.->> sig_prot_filter}

joined_df_name <- data.frame()

dataframe1 %>%
  filter(Gene %in% sig_prot_filter) %>%
  dplyr::select(Gene, Log2FC, FDR, Direction) %>%
  rbind(joined_df_name) %>%
  inner_join(dataframe2, by = c("Gene" = "Gene")) %>%
  dplyr::select(Gene, Log2FC.x, Log2FC.y, FDR.x, FDR.y, Direction.x, Direction.y) %>%
  return()
}

