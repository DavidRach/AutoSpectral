# check_critical.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.

#' @title check.critical
#'
#' @description Checks condition, if not true prints error message and stops execution.
#'
#' @param condition Essential condition(s) evaluated as a logical TRUE/FALSE.
#' @param error.msg Message to be returned to the user if condition is FALSE.
#' @return Returns error.msg and breaks if FALSE, no return if TRUE.
#' @export

check.critical <- function( condition, error.msg )
{
    if ( ! all( condition ) )
    {
        stop( sprintf( "%s\n", error.msg ), call. = FALSE )
    }
}

