#' JK Fourway Dataframe
#'
#' @param dataframe1 first dataframe, will correspond to "X" axis of four-way
#' @param dataframe2 second dataframe, will correspond to "Y" axis of four-way
#' @param gene_col column with gene names
#' @param fc_col column with fold change info to be plotted, usually Log2-FC or Log-FC
#' @param fdr_col column with p-value/q-value/FDR value info to be plotted
#' @param fc_thres threshold of fold change to consider. Default is 0 because preference is to still plot genes below the FC threshold as gray.
#' @param fdr_thres threshold of p/q/FDR value to plot, default 0.05
#' @return A data frame with shared significant genes/proteins across two contrasts, ready for plotting 4-ways
#' @export
#' @importFrom dplyr %>% count mutate case_when group_by inner_join arrange
#' @examples jk_sharedDE(testvol, fc_col = "log2FC", fdr_col = "qVal") #basic volcano plot using dataframe "testvol"
#' @examples jk_sharedDE(testvol, fc_col = "log2FC", fdr_col = "qVal", genenames = FALSE) #basic volcano plot using dataframe "testvol" with unlabeled points
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

#add column on direction based on above or below 0
dataframe1 = dataframe1 %>% dplyr::mutate(Direction = case_when(Log2FC > 0 ~ "UP", T ~ "DOWN"))
dataframe2 = dataframe2 %>% dplyr::mutate(Direction = case_when(Log2FC > 0 ~ "UP", T ~ "DOWN"))

dataframe1 %>%
  group_by(Gene) %>%
  count() %>% dplyr::filter(n>1) %>%
  dplyr::select(Gene) %>% unlist() %>%
  {.->>multiples.df1}

dataframe2 %>%
  group_by(Gene) %>%
  count() %>% dplyr::filter(n>1) %>%
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
  dplyr::filter(!Gene %in% multiples) %>%
  dplyr::filter(abs(Log2FC) > fc_thres,
         FDR < fdr_thres) %>%
  dplyr::select(Gene) %>% unlist() %>% {.->>sig_prot_a}


dataframe2 %>%
  dplyr::filter(!Gene %in% multiples) %>%
  dplyr::filter(abs(Log2FC) > fc_thres,
         FDR < fdr_thres) %>%
  dplyr::select(Gene) %>% unlist() %>% {.->>sig_prot_b}

unique(c(sig_prot_a, sig_prot_b)) %>% {.->> sig_prot_filter}

joined_df_name <- data.frame()

dataframe1 %>%
  dplyr::filter(Gene %in% sig_prot_filter) %>%
  dplyr::select(Gene, Log2FC, FDR, Direction) %>%
  rbind(joined_df_name) %>%
  inner_join(dataframe2, by = c("Gene" = "Gene")) %>%
  dplyr::select(Gene, Log2FC.x, Log2FC.y, FDR.x, FDR.y, Direction.x, Direction.y) %>%
  return()
}

#' JK Fourway Plot
#'
#' @param dataframe1 first dataframe, will correspond to "X" axis of four-way
#' @param dataframe2 second dataframe, will correspond to "Y" axis of four-way
#' @param fc_col column with fold change info to be plotted, usually Log2-FC or Log-FC
#' @param fdr_col column with p-value/q-value/FDR value info to be plotted
#' @param fc_thres threshold of fold change to consider. Default is 0 because preference is to still plot genes below the FC threshold as gray.
#' @param fdr_thres threshold of p/q/FDR value to plot, default 0.05
#' @param genenames boolean of whether gene names should be plotted on fourway
#' @param plot_fourway boolean of whether to return the fourway plot, or simply return the final dataframe generate
#' @param label_mode if labeling genes, what method (either [shared], showing shared significant genes, [unique], showing genes significant in only one contrast, [highlight], showing a select list of genes, or [all] genes. Default is "all".
#' @param color_mode affects if points are colored by either "overlap" (i.e., shared or unique between the two contrasts; default) or "quadrants" (includes direction per contrast as well)
#' @param symmetry boolean of whether the coordinates should resemble a square (i.e., same values above and below 0 per axis, obtained when param is TRUE) or not
#' @param color_axes boolean, default TRUE, of whether axes should be colored or not; uses the same colors used in color_mode quadrant
#' @param color_x_up for color_mode "quadrant", color for points defined as "up" in dataframe1
#' @param color_x_down for color_mode "quadrant", color for points defined as "down" in dataframe1
#' @param color_y_up for color_mode "quadrant", color for points defined as "up" in dataframe2
#' @param color_y_down for color_mode "quadrant", color for points defined as "down" in dataframe2
#' @param color_x for color_mode "overlap", color for points in dataframe 1
#' @param color_y for color_mode "overlap", color for points in dataframe 2
#' @param color_both_same for either color mode, sets color for genes significant in the same direction from both dataframes
#' @param color_both_diff for either color mode, sets color for genes significant in the different directions from both dataframes
#' @param color_ignore for either color mode, sets color for genes to be ignored (usually because they pass significance but are below FC thres)
#' @param max_overlaps setting for gene name labeling passed to geom_text_repel; higher value will allow for more genes to be plotted
#' @param contrast_1 name of the contrast in dataframe1, x axis label
#' @param contrast_2 name of the contrast in dataframe2, y axis label
#' @param label_size font size for gene name labels
#' @param highlight_list if label_mode is 'highlight', provide a list of gene/protein names to display
#' @param highlight_color if label_mode is 'highlight', what color to use for the points corresponding to gene/protein names provided in highlight_list
#' @param right_group name of experimental group on the right side (x > fc_thres) of x-axis
#' @param up_group  name of experimental group on the upper side (y > fc_thres) of y-axis
#' @param segment_trans transparency of the line that connects gene label to the plot of the gene, default is 0.3 (with 0 being fully transparent, and 1 fully opaque)
#' @param unity_alpha =1,
#' @param aspect aspect ratio of plot, default is 12/16,
#' @param point_size size of points when plotted, default 2
#' @param plot_title the title for the plot
#' @return A data frame with shared significant genes/proteins across two contrasts, ready for plotting 4-ways
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline element_text scale_y_continuous scale_x_continuous expansion annotate geom_abline xlab ylab xlim ylim theme element_blank scale_color_manual ggtitle
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% count mutate
#' @examples jk_4wayplot(df1, df2, fc_col = "log2FC", fdr_col = "qVal", genenames=FALSE) #basic fourway plot from two DE tables, no genes labeled

jk_4wayplot <- function(dataframe1,
                        dataframe2,
                        fc_col = "log2FC", #fold change column
                        fdr_col = "FDR", #p or q value column
                        fc_thres = 0.01,
                        fdr_thres = 0.05,
                        plot_fourway = TRUE,
                        genenames = TRUE,
                        label_mode = "all",
                        color_mode = "quadrants", #overlap or quadrants
                        symmetry = TRUE, #for axis symmetry
                        color_axes = TRUE,
                        color_x_up = "#4d528f",
                        color_x_down = "#a3a1c3",
                        color_y_up = "#de5b6c",
                        color_y_down = "#f3acb0",
                        color_x = "#9A9ABC",
                        color_y = "#F6A6A9",
                        color_both_same = "black",
                        color_both_diff = "blue",
                        color_ignore = "grey",
                        max_overlaps = 20,
                        contrast_1 = "Contrast1_Name",
                        contrast_2 = "Contrast2_Name",
                        label_size = 2.5, #gene name size
                        highlight_list = NULL,
                        highlight_color = "skyblue",
                        right_group = "group1",
                        up_group = "group2",
                        segment_trans = 0.3,
                        unity_alpha =1,
                        aspect = 12/16,
                        point_size = 2,
                        plot_title = "Significantly altered genes/proteins in both contrasts") {

  if (!color_mode %in% c("overlap", "quadrants")) stop("color_mode can only be either 'overlap' or 'quadrants'")

  fourway_df <- jk_sharedDE(dataframe1, dataframe2, fc_col = fc_col, fdr_col = fdr_col) %>%
    mutate(Sig.x = abs(Log2FC.x) > fc_thres & FDR.x < fdr_thres,
           Sig.y = abs(Log2FC.y) > fc_thres & FDR.y < fdr_thres) %>%
    mutate(color_overlap = case_when(Sig.x == TRUE & Sig.y == FALSE & abs(Log2FC.x) > fc_thres ~ "Sig_in_X",
                                     Sig.x == FALSE & Sig.y == TRUE & abs(Log2FC.y) > fc_thres ~ "Sig_in_Y",
                                     Sig.x == TRUE & Sig.y == TRUE & Log2FC.x * Log2FC.y > 0 ~ "SigBoth_same",
                                     Sig.x == TRUE & Sig.y == TRUE & Log2FC.x * Log2FC.y< 0 ~ "SigBoth_different",
                                     T ~ "Ignore")) %>%

    mutate(color_quadrants = case_when(Sig.x == TRUE & Sig.y == FALSE & Log2FC.x > fc_thres ~ "Up.X",
                                       Sig.x == TRUE & Sig.y == FALSE & Log2FC.x < -1*fc_thres ~ "Down.X",
                                       Sig.y == TRUE & Sig.x == FALSE & Log2FC.y > fc_thres ~ "Up.Y",
                                       Sig.y == TRUE & Sig.x == FALSE & Log2FC.y < -1*fc_thres ~ "Down.Y",
                                       Sig.x == TRUE & Sig.y == FALSE & Log2FC.x > fc_thres ~ "Up.X",
                                       Sig.x == TRUE & Sig.y == TRUE & Log2FC.x * Log2FC.y > 0 ~ "SigBoth_same",
                                       Sig.x == TRUE & Sig.y == TRUE & Log2FC.x * Log2FC.y < 0 ~ "SigBoth_different",
                                       T ~ "Ignore",
                                       abs(Log2FC.x) < fc_thres & abs(Log2FC.y) < fc_thres ~ "Ignore")
           ) %>%
    mutate(SigBothFilter = color_quadrants %in% c("SigBoth_same", "SigBoth_different")) %>%
    mutate(Label = color_quadrants %in% c("Up.X", "Down.X", "Up.Y", "Down.Y"))

  if(color_mode == "overlap"){
    fourway_df = fourway_df %>% arrange(color_overlap) #sort so that 'ignored' points plotted first
  }
  if(color_mode == "quadrants"){
    fourway_df = fourway_df %>% arrange(color_quadrants) #sort so that 'ignored' points plotted first
  }

  if (label_mode == "highlight") {
    fourway_df = fourway_df %>% dplyr::mutate(highlight_filter = Gene %in% highlight_list)
  }

  if(plot_fourway == FALSE){
    return(fourway_df)
    } else{

      four_way_plot_max.x <- fourway_df %>% dplyr::select(Log2FC.x) %>% max()
      four_way_plot_max.y <- fourway_df %>% dplyr::select(Log2FC.y) %>% max()
      four_way_plot_min.x <- fourway_df %>% dplyr::select(Log2FC.x) %>% min()
      four_way_plot_min.y <- fourway_df %>% dplyr::select(Log2FC.y) %>% min()
      n_sig_both_same <- fourway_df %>% dplyr::filter(color_overlap == "SigBoth_same") %>% nrow()
      n_sig_both_diff <- fourway_df %>% dplyr::filter(color_overlap == "SigBoth_different") %>% nrow()
      n_sig_x <- fourway_df %>% dplyr::filter(color_overlap == "Sig_in_X") %>% nrow()
      n_sig_y <- fourway_df %>% dplyr::filter(color_overlap == "Sig_in_Y") %>% nrow()
      n_sig_x_up <- fourway_df %>% dplyr::filter(color_quadrants == "Up.X") %>% nrow()
      n_sig_x_down <- fourway_df %>% dplyr::filter(color_quadrants == "Down.X") %>% nrow()
      n_sig_y_up <- fourway_df %>% dplyr::filter(color_quadrants == "Up.Y") %>% nrow()
      n_sig_y_down <- fourway_df %>% dplyr::filter(color_quadrants == "Down.Y") %>% nrow()

      four_way_plot_min <- min(four_way_plot_min.x, four_way_plot_min.y)
      four_way_plot_max <- max(four_way_plot_max.x, four_way_plot_max.y)

      xmargin.max =  four_way_plot_max.x + ((four_way_plot_max.x - four_way_plot_min.x)/15)
      xmargin.min = four_way_plot_min.x - ((four_way_plot_max.x - four_way_plot_min.x)/15)

      ymargin.max =  four_way_plot_max.y + ((four_way_plot_max.y - four_way_plot_min.y)/15)
      ymargin.min = four_way_plot_min.y - ((four_way_plot_max.x - four_way_plot_min.x)/15)

      #for symmetry, default
      abs_lim_x = max(abs(xmargin.min), abs(xmargin.max))
      abs_lim_y = max(abs(ymargin.min), abs(ymargin.max))

      max_lim_x = abs_lim_x
      min_lim_x = -1*abs_lim_x

      max_lim_y = abs_lim_y
      min_lim_y = -1*abs_lim_y

      if (symmetry ==FALSE){
        max_lim_x = xmargin.max
        min_lim_x = xmargin.min

        max_lim_y = ymargin.max
        min_lim_y = ymargin.min
      }

      plot_fourway <- fourway_df %>%
        ggplot(., aes(x=Log2FC.x, y=Log2FC.y, label=Gene)) +
        geom_hline(yintercept = (fc_thres), color="#898989", linetype = "dashed") +
        geom_hline(yintercept = -1*(fc_thres), color="#898989", linetype = "dashed") +
        geom_vline(xintercept = (fc_thres), color="#898989", linetype = "dashed") +
        geom_vline(xintercept = -1*(fc_thres), color="#898989", linetype = "dashed") +
        geom_abline(intercept=0, slope=1, color="#898989", alpha = unity_alpha) +
        xlim(c(four_way_plot_min,four_way_plot_max)) +
        ylim(c(four_way_plot_min,four_way_plot_max)) +
        ggplot2::theme_classic() +
        theme(legend.title = element_blank(),
              aspect.ratio = aspect,
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        scale_y_continuous(expand = expansion(mult = c(0,0),
                                              add = c(0,0)), #c(0,0),
                           limits = c(min_lim_y, max_lim_y)) +
        scale_x_continuous(expand = expansion(mult = c(0,0),
                                              add = c(0,0)), #c(0,0),
                           limits = c(min_lim_x, max_lim_x))

      if(color_mode == "overlap"){
        col_overlap =  c(color_ignore, color_x, color_y, color_both_same, color_both_diff)
        names(col_overlap) <- c("Ignore","Sig_in_X", "Sig_in_Y", "SigBoth_same", "SigBoth_different")

        plot_fourway = plot_fourway + geom_point(aes(color = color_overlap), size = point_size) +
          scale_color_manual(values = col_overlap)

        }
      if(color_mode=="quadrants") {
        col_quadrants = c(color_ignore, color_x_up, color_y_up, color_x_down, color_y_down, color_both_same, color_both_diff)
        names(col_quadrants) = c("Ignore", "Up.X", "Up.Y", "Down.X", "Down.Y", "SigBoth_same", "SigBoth_different")

        plot_fourway = plot_fourway + geom_point(aes(color = color_quadrants),
                                                 size = point_size) + scale_color_manual(values = col_quadrants)

      }

      if (color_axes==TRUE){
        plot_fourway = plot_fourway +

          #x-axis
          ggplot2::annotate(geom='segment', y = min_lim_y, yend = min_lim_y, x = (min_lim_x), xend = 0, linewidth=2, color = color_x_down) +
          ggplot2::annotate(geom='segment', y = min_lim_y, yend = min_lim_y, x = 0, xend = max_lim_x, linewidth=2, color = color_x_up) +

          #y-axis
          ggplot2::annotate(geom='segment', x = min_lim_x, xend = min_lim_x, y = (min_lim_y), yend = 0, linewidth=2, color = color_y_down) +
          ggplot2::annotate(geom='segment',x = min_lim_x, xend = min_lim_x, y = 0, yend = max_lim_y, linewidth=2, color = color_y_up)
        }


      if (genenames == TRUE){
        if (!label_mode %in% c("all", "shared", "unique", "highlight")) stop("label_mode can only be either 'all', 'shared', 'unique', or 'highlight'")

        if (label_mode == "highlight") {
          plot_fourway <- plot_fourway + ggplot2::geom_point(data=fourway_df[fourway_df$highlight_filter == T,],
                                                             aes(x = Log2FC.x, y = Log2FC.y),
                                                             color=highlight_color,
                                                             size = point_size, stroke =1)
          plot_fourway <- plot_fourway + ggrepel::geom_text_repel(aes(x = Log2FC.x, y = Log2FC.y,
                                                                      label = ifelse(highlight_filter == T, Gene,"")),
                                                                  size=label_size,
                                                                  max.overlaps = max_overlaps,
                                                                  segment.alpha = segment_trans,
                                                                  direction = "both",
                                                                  force = 1,
                                                                  min.segment.length = 0.25)
          }


        if(label_mode == "shared"){
          plot_fourway <- plot_fourway + geom_text_repel(data = fourway_df,
                                                         aes(x = Log2FC.x,
                                                             y = Log2FC.y,
                                                             fontface = "bold",
                                                             label = ifelse(SigBothFilter == TRUE , Gene,"")),
                                                         size=label_size,
                                                         max.overlaps = max_overlaps,
                                                         segment.alpha = segment_trans,
                                                         direction = "both",
                                                         force = 15,
                                                         min.segment.length = 0.25)
          }
        if (label_mode == "unique"){
        plot_fourway <- plot_fourway + geom_text_repel(data = fourway_df,
                                                           aes(x = Log2FC.x,
                                                               y = Log2FC.y,
                                                               fontface = "bold",
                                                               label = ifelse(Label == TRUE , Gene,"")),
                                                           size=label_size,
                                                           max.overlaps = max_overlaps,
                                                           segment.alpha = segment_trans,
                                                           direction = "both",
                                                           force = 15,
                                                           min.segment.length = 0.25)
        }

        if (label_mode == "all"){
          plot_fourway <- plot_fourway + geom_text_repel(data = fourway_df,
                                                         aes(x = Log2FC.x,
                                                             y = Log2FC.y,
                                                             fontface = "bold",
                                                             label = Gene),
                                                         size=label_size,
                                                         max.overlaps = max_overlaps,
                                                         segment.alpha = segment_trans,
                                                         direction = "both",
                                                         force = 15,
                                                         min.segment.length = 0.25)
        }

      }

      plot_fourway <- plot_fourway +
            xlab(paste0(contrast_1)) +
            ylab(paste0(contrast_2)) +
            ggtitle(paste0(plot_title))

      return(plot_fourway)
    }
  }

#unused code from prev versions

# labeling by color
# plot_fourway = plot_fourway + geom_point(aes(color = color_quadrants), size = point_size) +
#   scale_color_manual(values = c(color_x_up,
#                                 color_y_up,
#                                 color_both,
#                                 color_x_down,
#                                 color_y_down,
#                                 color_ignore),
#                      labels = c(paste0("Up in ", right_text, " (n=", n_sig_x_up, ")"),
#                                 paste0("Up in ", up_text," (n=",n_sig_y_up,")"),
#                                 paste0("Sig. in both (n=",n_sig_both,")"),
#                                 paste0("Down in ", right_text, " (n=", n_sig_x_down, ")"),
#                                 paste0("Down in ", up_text," (n=",n_sig_y_down,")"),
#                                 "")
#   ) +

# adding text
# +
#
#   #rightside
#   annotate(geom="text", x = xmargin.max,
#            label = paste0("Up in ", right_text),
#            y =four_way_plot_min.y*0.9,
#            size = 3,
#            hjust=1,
#            fontface =2,
#            color = color_x_up) +
#   #upper side
#   annotate(geom="text", y = ymargin.max,
#            label = paste0("Up in ", up_text),
#            x =four_way_plot_min.x*0.9,
#            size = 3,
#            hjust=0,
#            angle = 270,
#            fontface =2,
#            color = color_y_up)  +
#   #leftside
#   annotate(geom="text", x = four_way_plot_min.x,
#            label = paste0("  ",left_text),
#            y =four_way_plot_min.y*0.9,
#            size = 3,
#            hjust=0,
#            fontface =2,
#            color = color_x_down) +
#   #down side
#   annotate(geom="text", y = 0,
#            label = paste0(down_text, "  "),
#            x =four_way_plot_min.x*0.9,
#            size = 3,
#            hjust=0,
#            angle = 270,
#            fontface =2,
#            color = color_y_down)

# unused different way to plot points, couldnt figure how to make work with legends
# if (nrow(fourway_df %>% filter(color_overlap == "Ignore"))>0){
#   plot_fourway = plot_fourway +geom_point(data = subset(fourway_df, color_overlap == "Ignore"), color = color_ignore, size = point_size)
# }
# plot_fourway = plot_fourway +
# geom_point(data = subset(fourway_df, color_overlap == "Sig_in_X"), color = color_x, size = point_size) +
# geom_point(data = subset(fourway_df, color_overlap == "Sig_in_Y"), color = color_y, size = point_size) +
# geom_point(data = subset(fourway_df, color_overlap == "Sig_in_Both"), color = color_both, size = point_size) +
# scale_color_manual(values = c(color_x,
#                               color_y,
#                               color_both,
#                               color_ignore),
#                    labels = c(paste0("Sig. in ", contrast_1, " (n=", n_sig_x, ")"),
#                               paste0("Sig. in ", contrast_2," (n=",n_sig_y,")"),
#                               paste0("Sig. in both (n=",n_sig_both,")"),
#                               "")
#)

# if (nrow(fourway_df %>% filter(color_quadrants == "Ignore"))>0){
#   plot_fourway = plot_fourway +geom_point(data = subset(fourway_df, color_quadrants == "Ignore"), color = color_ignore, size = point_size)
# }
#plot_fourway = plot_fourway +
  # geom_point(data = subset(fourway_df, color_quadrants == "Up.X"), color = color_x_up, size = point_size) +
  # geom_point(data = subset(fourway_df, color_quadrants == "Up.Y"), color = color_y_up, size = point_size) +
  # geom_point(data = subset(fourway_df, color_quadrants == "Down.X"), color = color_x_down, size = point_size) +
  # geom_point(data = subset(fourway_df, color_quadrants == "Down.Y"), color = color_y_down, size = point_size) +
  # geom_point(data = subset(fourway_df, color_quadrants == "Sig_in_Both"), color = color_both, size = point_size) +

# 2025-8-5 new unused code

# mutate(color_overlap = factor(color_overlap, levels=c("Ignore",
#                                                       "Sig_in_X",
#                                                       "Sig_in_Y",
#                                                       "Sig_in_Both"
#                                                       ))
#        ) %>%
# mutate(color_quadrants = factor(color_quadrants, levels = c("Ignore",
#                                                             "Up.X",
#                                                             "Up.Y",
#                                                             "Down.X",
#                                                             "Down.Y",
#                                                             "Sig_in_Both"
#                                                             ))
#        ) %>%

# plot_fourway = plot_fourway +
#   geom_point(aes(color=color_quadrants),size = point_size) +
#   scale_color_manual(values = c(color_ignore,
#                                 color_x_up,
#                                 color_y_up,
#                                 color_x_down,
#                                 color_y_down,
#                                 color_both
#                                 ),
#                      labels = c("Below FC/FDR thresholds",
#                                 paste0("Up in ", right_group, " (n=", n_sig_x_up, ")"),
#                                 paste0("Up in ", up_group," (n=",n_sig_y_up,")"),
#                                 paste0("Down in ", right_group, " (n=", n_sig_x_down, ")"),
#                                 paste0("Down in ", up_group," (n=",n_sig_y_down,")"),
#                                 paste0("Sig. in both (n=",n_sig_both,")")
#                                 )
#                      )


# plot_fourway = plot_fourway +
#   geom_point(aes(color=color_overlap),size = point_size) +
#   scale_color_manual(values = c(color_ignore,
#                                 color_x,
#                                 color_y,
#                                 color_both
#                                 ),
#                      labels = c("Below FC/FDR thresholds",
#                                 paste0("Sig. in ", contrast_1, " (n=", n_sig_x, ")"),
#                                 paste0("Sig. in ", contrast_2," (n=",n_sig_y,")"),
#                                 paste0("Sig. in both (n=",n_sig_both,")")
#                                 )
#                      )
