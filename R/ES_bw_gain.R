#' @title ES Weight Gain Calculator
#' @description This function can be used to calculate pup weight gain in early-life stress (ES) experiments, by subtracting average postnatal day (P)2 weights from the final P9 weights, treating male and female pups separately.
#'
#' @param df Dataframe with body weights
#' @param sex Name of column with pup sex info, default is "Sex"
#' @param bwp2 Name of column with P2 weights, default is "BWP2"
#' @param bwp9 Name of column with P9 weights, default is "BWP9"
#' @param nest Name of column with dam/nest info, default is "Dam"
#' @param pupID Name of column with pup IDs, default is "PupID"
#' @param condition Name of column with condition info, default is "Condition"
#' @param delta_name Name to use for column with body weight gain, default is "BW_gain"
#' @param analysis_mode Whether to perform this on a pup or nest level
#'
#' @return Returns a dataframe of P9 - P2 body weights either at the pup or nest level.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%#'
#' @examples #generate test column with simulated values
#' @examples test_df <- data.frame(PupID = rep("Test",10),Condition = rep("CTL", 10), BWP2 = runif(10, 1.5, 3.5), BWP9 = runif(10, 4, 5.2), Sex = c(rep("m", 5),  rep("f", 5)), Dam = c(rep("F01",3), rep("F02",3), rep("F01",2), rep("F02",2)))
#' @examples ES_bw_gain(test_df) #generates a column BW_gain of P9 - avgP2 weights per pup
#'
#'
ES_bw_gain <- function(df, #dataframe
                       sex = "Sex", #col defining sex of pups
                       bwp2 = "BWP2", #col with body weights at P2
                       bwp9 = "BWP9", #col with body weights at P9
                       nest = "Dam", #col with nest information
                       pupID = "PupID", #col with pup IDs
                       condition = "Condition", #col with condition info
                       delta_name = "BW_gain", #name to use for column containing result of body weight subtraction
                       analysis_mode = c('pups', 'nest')

){

  analysis_mode <- match.arg(analysis_mode)

  #end solution: rename columns to fit what is in the code
  names(df)[names(df) == {{bwp2}}] <- "BWP2"
  names(df)[names(df) == {{bwp9}}] <- "BWP9"
  names(df)[names(df) == {{sex}}] <- "Sex"
  names(df)[names(df) == {{nest}}] <- "Dam"
  names(df)[names(df) == {{condition}}] <- "Condition"
  names(df)[names(df) == {{pupID}}] <- "PupID"

  if (analysis_mode == 'pups'){
    calculated <- df %>%
      dplyr::group_by(Dam, Sex) %>%

      dplyr::summarize(PupID = PupID,
                       Condition = Condition,
                       Dam = Dam,
                       BWP2 = BWP2,
                       P2_Avg = mean(BWP2, na.rm = TRUE),
                       BWP9 = BWP9,
                       bw_delta = BWP9 - mean(BWP2, na.rm = TRUE))
  }

  if (analysis_mode == 'nest'){
    calculated <- df %>%
      dplyr::group_by(Dam, Sex) %>%

      dplyr::summarize(Dam = Dam,
                       Condition = Condition,
                       avgBWP2 = mean(BWP2, na.rm = TRUE),
                       avgBWP9 = mean(BWP9, na.rm = TRUE),
                       bw_delta = avgBWP9 - avgBWP2)
  }

  names(calculated)[names(calculated) == "bw_delta"] <- {{delta_name}}

  return(calculated)

}


