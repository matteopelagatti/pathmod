#' Russett Dataset
#'
#' Data set from Russett (1964) about agricultural inequality, industrial development
#' and political instability.
#'
#' @format A data frame with 47 rows and 11 variables. The variables may be used to construct three latent concepts: 1) AGRIN=Agricultural Inequality, 2) INDEV=Industrial Development, 3) POLINS=Political Instability.
#' \describe{
#'   \item{gini}{Inequality of land distribution (AGRIN)}
#'   \item{farm}{Percentage of farmers that own half of the land (AGRIN)}
#'   \item{rent}{Percentage of farmers that rent all their land}
#'   \item{gnpr}{Gross national product per capita (INDEV)}
#'   \item{labo}{Percentage of labor force employed in agriculture (INDEV)}
#'   \item{inst}{Instability of executive (1945-1961) (POLINS)}
#'   \item{ecks}{Number of violent internal war incidents (1941-1961) (POLINS)}
#'   \item{death}{Number of people killed as a result of civic group violence (1950-1962) (POLINS)}
#'   \item{demostab}{Political regime: stable democracy (POLINS)}
#'   \item{demoinst}{Political regime: unstable democracy (POLINS)}
#'   \item{dictator}{Political regime: dictatorship (POLINS)}
#' }
#' @references Russett B.M. (1964) Inequality and Instability: The Relation of Land Tenure to Politics. \emph{World Politics} 16:3, 442-454.
#' @examples
#' data(russett)
#' russett
"russett"
