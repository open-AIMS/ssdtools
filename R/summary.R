# Copyright 2015-2023 Province of British Columbia
# Copyright 2021 Environment and Climate Change Canada
# Copyright 2023-2024 Australian Government Department of Climate Change,
# Energy, the Environment and Water
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' @export
summary.tmbfit <- function(object, ...) {
  chk_unused(...)
  x <- list()
  x$dist <- .dist_tmbfit(object)
  x$estimates <- estimates(object)
  class(x) <- "summary_tmbfit"
  x
}

#' @export
summary.fitdists <- function(object, ...) {
  chk_unused(...)
  x <- list()
  x$fits <- lapply(object, summary)
  x$censoring <- .censoring_fitdists(object)
  x$nrow <- nrow(.data_fitdists(object))
  x$rescaled <- .rescale_fitdists(object)
  x$weighted <- .weighted_fitdists(object)
  x$unequal <- .unequal_fitdists(object)
  x$min_pmix <- .min_pmix_fitdists(object)
  class(x) <- "summary_fitdists"
  x
}
