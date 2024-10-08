% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scales.R
\name{ssd_label_comma}
\alias{ssd_label_comma}
\title{Label numbers with significant digits and comma}
\usage{
ssd_label_comma(digits = 3, ..., big.mark = ",")
}
\arguments{
\item{digits}{A whole number specifying the number of significant figures.}

\item{...}{Unused.}

\item{big.mark}{A string specifying used between every 3 digits to separate thousands on the x-axis.}
}
\value{
A "labelling" function that takes a vector x and
returns a character vector of \code{length(x)} giving a label for each input value.
}
\description{
Label numbers with significant digits and comma
}
\examples{
ggplot2::ggplot(data = ssddata::anon_e, ggplot2::aes(x = Conc / 10)) +
  geom_ssdpoint() +
  ggplot2::scale_x_log10(labels = ssd_label_comma())
}
\seealso{
\code{\link[scales:label_number]{scales::label_comma()}}
}
