#' Replication dataset from Del Ponte (2020)
#'
#' Ponte, Alessandro Del.
#' "The influence of foreign elite rhetoric:
#' National identity, emotions, and attitudes toward austerity."
#' European Union Politics 22.1 (2021): 155-178.
#'
#' A dataset containing the results of an online survey
#' experiment conducted on Italian respondents who read an article about 
#'
#' @format A data frame with 707 rows and 53 variables: \describe{
#'   \item{angry_bin}{Self-reported anger after reading the the article}
#'   \item{itaid_bin}{Binary measure of strength of Italian identity}
#'   \item{t_commonality}{t=Treatment indicator for if the article praises (1)
#'   or blames (0) Italy} \item{north}{Indicator for respondent living in
#'   Northern Italy} \item{satisf}{Answer on a 4-point scale (rescaled to 0-1)
#'   to the question
#'   "In general, how satisfied or dissatisfied are you with your economic situation?"}
#'   \item{sopscale}{Measure of how politically sophisticated the respondent is
#'   by correct answers to 2 factual question. Rescaled to 0-1.}
#'   \item{Corriere}{Indicator for if the respondent is a regular reader of the
#'   Corriere blog where the survey was advertised.}}
#' @source \url{https://journals.sagepub.com/doi/suppl/10.1177/1465116520966653}
"delponte"
