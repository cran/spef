\name{skiTum}
\docType{data}
\alias{skiTum}
\title{Skin cancer chemoprevention trial}
\usage{skiTum}
\description{ A data frame contains data on the recurrence of two types of
  non-melanoma skin tumor, basal cell carcinomas and squamous cell carcinomas.
  The data was obtained from a double-blinded, placebo-controlled,
  randomized Phase III clinical trial (Bailey et al., 2010) conducted at
  the University of Wisconsin Comprehensive Cancer Center.
  This data set contains records of 61 elder female patients (60 years or older)
  with a history of skin cancer randomized into two groups:
  a treatment group with oral difluoromethylornithine (DFMO) at a daily dose of 0.5 gram/m2 and a placebo group.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{id}}{: patient id (repeated for each recurrence)}
    \item{\code{time}}{: observation time}
    \item{\code{age}}{: patient's age at enrollment}
    \item{\code{male}}{: gender; male = 1, female = 0}
    \item{\code{dfmo}}{: treatment (DFMO) group = 1; placebo = 0}
    \item{\code{priorTumor}}{: number of prior tumor from diagnosis to randomization}
    \item{\code{countBC}}{: number of new basal cell carcinomas tumors since last observation time}
    \item{\code{countSC}}{: number of new squamous cell carcinomas
      tumors since last observation time}
    \item{\code{count}}{: number of new tumors since last observation
      time (countBC + countSC)}
  }
}
\references{
  Bailey, H. H., Kim, K., Verma, A. K., Sielaff, K., Larson, P. O.,
  Snow, S., Lenaghan, T., Viner, J. L., Douglas, J., and  Dreckschmidt,
  N. E. (2010). A randomized, double-blind, placebo-controlled phase 3
  skin cancer prevention study of alpha-difluoromethylornithine in
  subjects with previous history of skin cancer. \emph{Cancer Prevention
    Research}, \bold{3}, 35--47.
}
\seealso{
  \code{\link{skinTumor}}
}
\examples{
data(skiTum)
library(ggplot2)
data(skiTum)
skiTum$dfmo <- factor(skiTum$dfmo, levels = c(1, 0), labels = c("placebo", "DFMO"))
ggplot(skiTum, aes(time, id, height = 2, width = 15)) +
    geom_tile(aes(fill = count)) + theme_bw() +
    facet_grid(dfmo ~ ., scales = "free_y") + 
    scale_fill_gradient(low = "grey", high = "black") + labs(fill="Count") +
    scale_x_continuous(breaks = seq(0, 1900, 200)) + xlab("Time in days")
}
\keyword{survival}
\keyword{panel count}
