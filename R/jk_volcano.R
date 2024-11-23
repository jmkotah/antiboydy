#' JK Volcano Plot v1
#'
#' @param df dataframe to be plotted
#' @param fc_col column with fold change info to be plotted, usually Log2-FC or Log-FC
#' @param fdr_col column with p-value/q-value/FDR value info to be plotted
#' @param fc_thres threshold of fold change to plot
#' @param fdr_thres threshold of p/q/FDR value to plot, default 0.05
#' @param contrast_name optional: name of contrast (used in title)
#' @param FC_log log base of the foldchange column, used in the x-axis plot (default is 2)
#' @param aspect aspect ratio to use, can be expressed as decimal or fraction; default is 1
#' @param color_left color of dots on the left side of the volcano plot
#' @param color_right color of dots on the right side of the volcano plot
#' @param point_size size of dots in plot
#' @param label_size size of labels in plot
#' @param xmargin added to the color
#' @param genenames TRUE/FALSE, whether gene/protein names should be displayed next to dots, default is TRUE
#' @param max_overlap Maximum number of gene names to try to plot, based on geom_text_repel
#' @param segment_trans Transparency of line that connects gene names to dots, default 0.3, based on geom_text_repel
#' @param right_text Show text on lower right part of plot, uses "color_right"
#' @param left_text Show text on lower left part of plot, uses "color_left"
#' @param contrast_text_size size of text from 'right_text' and 'left_text'
#' @return A volcano plot showing significantly up/downregulated genes or proteins
#' @export
#' @importFrom ggplot2 aes geom_point geom_hline geom_vline element_text scale_y_continuous scale_x_continuous expansion geom_segment annotate xlab
#' @importFrom ggrepel geom_text_repel
#' @examples jk_volcano(testvol, fc_col = "log2FC", fdr_col = "qVal") #basic volcano plot using dataframe "testvol"
#' @examples jk_volcano(testvol, fc_col = "log2FC", fdr_col = "qVal", genenames = FALSE) #basic volcano plot with unlabeled points
jk_volcano <- function(df, #dataframe
                       fc_col, #fold change column
                       fdr_col, #p or q value column
                       fc_thres = 0.1, #foldchange cutoff
                       fdr_thres = 0.05, #fdr threshold for significance
                       #filter #column with label to determine filtering,
                       contrast_name = "Insert Contrast Here",
                       FC_log = 2, #base-log foldchange is expressed in,
                       aspect = 1, #aspect ratio to use
                       color_left = "#de5b6c",
                       color_right = "#4d528f",
                       point_size = 1.25, #size of volcano points
                       label_size = 2.5,
                       xmargin = 0.1,
                       #ymargin = 0, unused in current version
                       genenames = TRUE,
                       max_overlap = 10,
                       segment_trans = 0.3,
                       right_text = "Contrast_Rt",
                       left_text = "Contrast_Lt",
                       contrast_text_size =4
) {

  #code developed for Log2FC and FDR, so renaming everything as such; in principle works with Log10FC, or p-val instead of FDR
  names(df)[names(df) == {{fc_col}}] <- "Log2FC"
  names(df)[names(df) == {{fdr_col}}] <- "FDR"

  df %>%
    dplyr::filter(FDR<0.05 & abs(Log2FC) > fc_thres) %>%
    dplyr::select(Log2FC) %>% abs() %>% max() %>% 
    {.->> volcano_plot_limits}

  df %>%
    dplyr::filter(FDR<0.05 & abs(Log2FC) > fc_thres) %>%
    dplyr::select(FDR) %>% min() %>% -log10(.) %>%
    {.->> pval_limit}

  df %>%
    dplyr::mutate(filter = FDR<fdr_thres & abs(Log2FC) > fc_thres) %>%
    dplyr::mutate(filter1 = FDR<fdr_thres & Log2FC < -1* fc_thres) %>% #negative filter
    dplyr::mutate(filter2 = FDR<fdr_thres & Log2FC > 1* fc_thres) %>% #positive filter
    {.->>test.df}

  volcano_plot <- ggplot2::ggplot(test.df) +
    geom_point(data = subset(test.df, filter == FALSE),
                        aes(x=Log2FC, y=-log10(FDR)), color = "gray", size = point_size) +
    geom_point(data = subset(test.df, filter1 == TRUE),
                        aes(x=Log2FC, y=-log10(FDR)), color = color_left, size = point_size) +
    geom_point(data = subset(test.df, filter2 == TRUE),
                        aes(x=Log2FC, y=-log10(FDR)), color = color_right, size = point_size) +

    geom_hline(yintercept=-log10(fdr_thres), linetype="dashed", color = "grey") +
    geom_vline(xintercept=(fc_thres), linetype="dashed", color = "grey") +
    geom_vline(xintercept=-1*(fc_thres), linetype="dashed", color = "grey")



  if (genenames == TRUE){
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(aes(x = Log2FC,
                          y = -log10(FDR),
                          label = ifelse(filter1 == T, Gene,"")),
                      size=label_size,
                      max.overlaps = max_overlap,
                      segment.alpha = segment_trans,
                      #nudge_x = -0.25,
                      direction = "both",
                      force = 1,
                      xlim = c(NA,fc_thres/-2),
                      ylim = c(1.3, NA),
                      #force_pull = 0,
                      min.segment.length = 0.25) +

      ggrepel::geom_text_repel(aes(x = Log2FC,
                          y = -log10(FDR),
                          label =  ifelse(filter2 == T, Gene,"")),
                      size=label_size,
                      max.overlaps = max_overlap,
                      segment.alpha = segment_trans,
                      #nudge_x = 0.25,
                      direction = "both",
                      force = 1,
                      xlim = c(fc_thres/2, NA),
                      ylim = c(1.3, NA),
                      #force_pull = 1,
                      min.segment.length = 0.25) }


  volcano_plot <- volcano_plot +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          aspect.ratio = aspect) +

    ggplot2::xlab(expression(paste(Log[2], " Fold Change"))) +
    ggplot2::ylab(expression(paste("-",Log[10]," (FDR corrected p-value)"))) +
    ggplot2::xlim(c(-1*(volcano_plot_limits+xmargin),volcano_plot_limits+xmargin)) +
    ggplot2::ggtitle(paste0(contrast_name))

  if (FC_log == 10){
    volcano_plot <- volcano_plot +
      xlab(expression(paste(Log[10], " Fold Change")))
  }


  volcano_plot <- volcano_plot +
    scale_y_continuous(expand = expansion(mult = c(0,0),
                                          add = c(0,0)), #c(0,0),
                       limits = c(0,1.05*pval_limit)) +
    scale_x_continuous(expand = expansion(mult = c(0,0),
                                          add = c(0,0))
    ) +

    geom_segment(aes(y = 0, yend = 0, x = -(volcano_plot_limits+ xmargin), xend = 0), size=2, color = color_left) +
    geom_segment(aes(y = 0, yend = 0, x = 0, xend = volcano_plot_limits + xmargin), size=2, color = color_right) +

    #rightside
    annotate(geom="text", x = volcano_plot_limits + xmargin/2,
             label = right_text,
             y =pval_limit*.03,
             size = contrast_text_size,
             hjust=1,
             fontface =2,
             color = color_right) +

    #leftside
    annotate(geom="text", x = -1*(volcano_plot_limits + xmargin/2),
             label = left_text,
             y =pval_limit*.03,
             size = contrast_text_size,
             hjust=0,
             fontface =2,
             color = color_left)

  volcano_plot
}
