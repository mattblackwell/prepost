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
#' @source <https://journals.sagepub.com/doi/suppl/10.1177/1465116520966653>
"delponte"


#' Replication dataset from Horowitz and Klaus (2020)
#'
#' Horowitz, J., Klaus, K. "Can Politicians Exploit Ethnic Grievances?
#' An Experimental Study of Land Appeals in Kenya." Political Behavior
#' 42, 35–58 (2020). https://doi.org/10.1007/s11109-018-9485-1
#'
#' A dataset containing the results of an experimental study of political
#' appeals about land greivances on candidate support in Kenya
#'
#' @format A data frame with 375 rows and 9 variables: \describe{
#'   \item{support}{Support for the hypothetical candidate making the appeal}
#'   \item{land_insecure}{Moderator of interest measuring if the respondent
#' is land insecure}
#'   \item{treat_comb}{Treatment indicator for land-based appeals}
#'   \item{prepost}{Indicator for if the moderator was measured before (0) or after(1) treatment}
#' \item{age}{Age of the respondent}
#'   \item{female}{Indicator for if the respondent identifies as female}
#'   \item{close_own}{Indicator for if the respondent feels close to their own ethnic group}
#'   \item{educ}{Level of education (1=no schooling, 2 = some primary,
#'     3 = complete primary, 4 = some secondary, 5 = completed secondary,
#'     6 = college, 7 = some university, 8 = completed university)},
#'   \item{treat}{Original treatment variable indicating control (0),land-based appeal (1),
#' and land- and ethnic-based appeals (2).}
#' }
#' @source \doi{10.7910/DVN/XBWR8N}
"land_experiment"

