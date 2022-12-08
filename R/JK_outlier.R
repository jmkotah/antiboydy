#' @title IQR-based Outlier Detection (JMK)
#' @description A function to show outliers from a set of numbers as tested using the IQR method, either by displaying the outlier values when present, or returning the upper and lower criteria for deciding if something is an outlier.
#'
#' @param variable A list of numeric values to subject to outlier testing
#' @param print_mode Default option "values" will return actual outliers if present; "cutoffs" with  return criteria to determine outliers
#' @param multiplier Number multiplied to the IQR to set outlier thresholds; default is 1.5
#' @param cutoff_type For use when print mode set to 'cutoffs'. Choose to display the higher ('high') or lower ('low') cutoffs, or "both"
#'
#' @return Returns either a list of values that are outliers (and their relative position in the query), or numbers to use to filter data for outliers
#' @export
#'
#' @examples #returns 15 and -10 as outliers
#' @examples outlier.jk(c(1,2,3,1,3,4,15,-10))
#'
#' @examples #shows the range values must be in to be deemed 'not an outlier'
#' @examples outlier.jk(c(1,2,3,1,3,4,15,-10) , print_mode = "cutoffs", cutoff_type = "both")

outlier.jk <- function(variable,
                       print_mode = c('values', 'cutoffs'),
                       multiplier = 1.5,
                       #print_outliers = TRUE,
                       cutoff_type = c("high", "low", "both")
                       #print_high = FALSE,
                       #print_low = FALSE
                       ) {
  stopifnot(is.numeric(variable))
  #stopifnot(print_mode %in% c('values', 'cutoffs'))
  print_mode <- match.arg(print_mode)
  cutoff_type <- match.arg(cutoff_type)

    Q1<-summary(variable)[[2]]
    Q3<-summary(variable)[[5]]

    mult <- (Q3-Q1)*multiplier
    low <- Q1-mult
    high <- Q3+mult

    index<-c(which(variable<low),which(variable>high))
    index<-sort(index, decreasing=FALSE)
    #print(variable[index])
    #return(length(index)>0)


    if (print_mode == "values") {
      if (length(index)>0){
        x<-paste0("Outlier: ", variable[index], " in position ", index, " (", multiplier,"*IQR Method)", collapse= '\n')
        cat(x)
      } else{
        writeLines("No outliers. (1.5IQR Method)")
      }
    }

    if (print_mode == "cutoffs") {
      if (cutoff_type == "high"){
        cat(paste0("Higher cutoff is: ", high))
      }

      if (cutoff_type == "low"){
        cat(cat(paste0("Lower cutoff is: ", low)))
      }
      if (cutoff_type == "both"){
        cat(paste0("Higher cutoff is: ", high, " and lower cutoff is: ", low))
      }
    }

  }
