# fit_robust_linear_model.r

#' Fit Robust Linear Model
#'
#' This function returns a matrix by rows, with the intercept and p-value,
#'     and coefficient and p-value, of a robust linear model fitted to the
#'     input data. It reverts to a standard linear model in case of no convergence.
#'
#' @title Fit Robust Linear Model
#' @description Returns a matrix by rows, with the intercept and p-value, and
#'     coefficient and p-value, of a robust linear model fitted to the input data.
#'     Reverts to a standard linear model in case of no convergence.
#'
#' @importFrom MASS rlm
#' @importFrom stats lm
#' @param x.data A vector containing the predictor variable data.
#' @param y.data A vector containing the response variable data.
#' @param x.name The name of the predictor variable.
#' @param y.name The name of the response variable.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'     A list of essential parameters.
#' @param fix.unmix Logical, default is FALSE. If TRUE, sets coefficient to zero
#'     in case of failed convergence. Used for fix.my.unmix.
#' @return A matrix with the intercept and p-value, and coefficient and p-value.
#' @export



fit.robust.linear.model <-  function( x.data, y.data, x.name, y.name, asp,
                                      fix.unmix = FALSE )
{
  X.data.int <- cbind( 1, x.data )
  xy.model <- rlm( X.data.int, y.data, maxit = asp$rlm.iter.max )

    if ( xy.model$converged ) {

      xy.coef <- xy.model$coefficients

    } else if ( ! xy.model$converged & fix.unmix ) {

      xy.coef <- 0

    } else {
      cat( sprintf( "WARNING: rlm of %s ~ %s did not converge - using ols instead\n",
                    y.name, x.name ), file = stderr() )

      xy.data <- data.frame( x = x.data, y = y.data )

      xy.model <- lm( y ~ x, xy.data )

      xy.coef <- xy.model$coefficients
    }

  dimnames( xy.coef ) <- NULL
  xy.coef
}


