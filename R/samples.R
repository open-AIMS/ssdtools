xsample_estimates <- function(x, args, what) {
  if (grepl("^ssd_p", what)) {
    args$q <- x
  } else {
    args$p <- x
  }
  samples <- do.call(what, args)
  samples
}

sample_estimates <- function(estimates, what, x, .names = NULL) {
  args <- transpose(estimates, .names = .names)
  args <- purrr::map_depth(args, 2, function(x) {
    if (is.null(x)) NA_real_ else x
  })
  args <- lapply(args, as.double)
  x <- lapply(x, xsample_estimates, args, what)
  x
}
