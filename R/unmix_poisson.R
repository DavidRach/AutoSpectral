# unmix_poisson.r

#' @title Unmix Using Poisson Regression
#'
#' @description This function performs unmixing of raw data using Poisson regression,
#'     with iterative reweighted least squares (IRLS) and fallback methods for
#'     cells that fail to converge.
#'
#' @importFrom stats glm glm.control poisson coef
#' @importFrom MASS rlm
#'
#' @param raw.data Matrix containing raw data to be unmixed.
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param allow.negative Logical indicating whether to allow negative coefficients.
#'
#' @return A matrix containing the unmixed data.
#' @export


unmix.poisson <- function( raw.data, spectra, asp, allow.negative = TRUE ) {

  # Initialize with OLS unmixing
  X <- t( spectra )
  coef.matrix <- unmix.ols( raw.data, spectra )

  # Handle zeros
  raw.data[ raw.data <= 0] <- 1e-6
  coef.matrix[ coef.matrix <= 0] <- 1e-6

  # Convergence tolerance
  tol <- 1e-6
  convergence <- rep( FALSE, nrow( raw.data ) )

  # Track cells that fail to converge for fallback method
  problem.cells <- integer( 0 )

  # IRLS implementation
  for ( iter in 1:asp$rlm.iter.max ) {
    # Skip cells that have already converged
    if ( all( convergence ) ) break

    # Store previous coefficients to check convergence
    prev.coef <- coef.matrix

    # Calculate predicted values (mu) for all cells simultaneously
    mu <- coef.matrix %*% t( X )  # cells × detectors

    # Ensure mu is positive
    mu[ mu <= 0 ] <- 1e-6

    # Process in batches for memory efficiency
    batch.size <- min( 5000, nrow( raw.data ) )
    batch.n <- ceiling( nrow( raw.data )/ batch.size )

    for ( b in 1 : batch.n ) {
      start.idx <- ( b-1 ) * batch.size + 1
      end.idx <- min( b*batch.size, nrow( raw.data ) )
      batch.idx <- start.idx : end.idx

      # Skip batch if all cells have converged
      if ( all( convergence[ batch.idx ] ) ) next

      # Non-converged cells in this batch
      active.idx <- batch.idx[ !convergence[ batch.idx ] ]

      # Process each non-converged cell in batch
      for ( j in 1 : length( active.idx ) ) {
        i <- active.idx[ j ]
        y_i <- raw.data[ i, ]
        mu_i <- mu[ i, ]

        # Calculate weights
        w_i <- 1 / mu_i
        w_i[ !is.finite( w_i ) ] <- 1e6  # Handle very small mu values

        # Weighted design matrix
        XtW <- t( X * w_i )

        # Form normal equations
        XtWX <- XtW %*% X
        XtWz <- XtW %*% y_i

        # Solve using QR decomposition
        tryCatch({
          new.coef <- qr.solve( XtWX, XtWz )

          # Handle negative coefficients based on allow.negative parameter
          if ( !allow.negative ) {
            new.coef[ new.coef < 0 ] <- 1e-6
          }

          # Update coefficients
          coef.matrix[ i, ] <- new.coef
        },
        error = function( e ) {
          # Mark this cell for fallback method
          problem.cells <- c( problem.cells, i )
        })
      }
    }

    # Check convergence for all cells
    diffs <- rowSums( ( coef.matrix - prev.coef )^2 ) / rowSums( prev.coef^2 + 1e-10 )
    convergence <- diffs < asp$rs.delta.threshold.change

    # Break if overall change is small
    if ( max( diffs[ !convergence ], 0) < asp$rs.delta.threshold.change * 10 ) break
  }

  # Handle problem cells with fallback methods
  problem.cells <- unique( c( problem.cells, which( !convergence ) ) )

  if ( length( problem.cells ) > 0) {
    for ( i in problem.cells ) {

      y_i <- raw.data[ i, ]
      coef.matrix_i <- coef.matrix[ i, ]
      coef.matrix_i[ coef.matrix_i <= 0] <- 1e-6

      # Try glm first
      fit <- try(
        suppressWarnings(
          glm( y_i ~ X - 1,
              family = poisson( link = "identity" ),
              start = coef.matrix_i,
              control = glm.control( maxit = asp$rlm.iter.max ),
              silent = TRUE )
        )
      )

      if ( !inherits(fit, "try-error") && fit$converged ) {
        coef.matrix[ i, ] <- coef( fit )

        # Handle negative coefficients based on allow.negative parameter
        if ( !allow.negative ) {
          coef.matrix[ i, coef.matrix[ i, ] < 0] <- 1e-6
        }
      } else {
        # Try rlm as fallback
        fit <- try(
          MASS::rlm( y_i ~ X - 1,
                    init = coef.matrix[ i, ],
                    maxit = asp$rlm.iter.max )
        )

        if ( !inherits( fit, "try-error" ) ) {
          coef.matrix[ i, ] <- coef( fit )

          # Handle negative coefficients based on allow.negative parameter
          if ( !allow.negative ) {
            coef.matrix[ i, coef.matrix[ i, ] < 0 ] <- 1e-6
          }
        }
        # If both fail, keep the last iteration's values
      }
    }
  }

  colnames( coef.matrix ) <- rownames( spectra )

  return( coef.matrix )
}
