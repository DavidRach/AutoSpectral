# check_critical.r

#' @title check.critical
#'
#' @description Checks condition, if not true prints error message and
#' stops execution.
#'
#' @param condition Essential condition(s) evaluated as a logical TRUE/FALSE.
#' @param error.msg Message to be returned to the user if condition is FALSE.
#'
#' @return Returns error.msg and breaks if FALSE, no return if TRUE.
#' @export

check.critical <- function( condition, error.msg )
{
    if ( ! all( condition ) )
    {
        stop( sprintf( "%s\n", error.msg ), call. = FALSE )
    }
}

