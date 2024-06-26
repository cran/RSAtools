% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportRSA.bootstrap.R
\name{exportRSA.bootstrap}
\alias{exportRSA.bootstrap}
\title{Export bootstrapped parameters of polynomial model}
\usage{
exportRSA.bootstrap(RSAbootstrap_object)
}
\arguments{
\item{RSAbootstrap_object}{A matrix output generated by lavaan: \code{bootstrapLavaan(RSA_object$models$name_final, FUN="coef")}}
}
\value{
A table of mean values and 95% confidence intervals of polynomial parameters
}
\description{
Utilitary function to export bootstrapped parameters, as part of STEP3 of the 3-step identification strategy (see \code{RSAmodel}).
Utilitary function to export bootstrapped parameters
}
\examples{
###### Export 95\% CI of bootstrapped estimates of polynomial
#Estimate a model: FM26_PARALLELASYMWEAK (simulation data)
RSA_NSfit  <- RSAmodel(formula= engagement ~ needs*supplies,
data= sim_NSfit, model= c("FM26_PARALLELASYMWEAK"))
#Bootstrapped sampling with lavaan
RSA_NSfit_boot  <- lavaan::bootstrapLavaan(RSA_NSfit$models$FM26_PARALLELASYMWEAK,
 R= 10,FUN="coef")
#Export results in a table
RSA_NSfit_boot_exp <- exportRSA.bootstrap(RSA_NSfit_boot)
RSA_NSfit_boot_exp
}
