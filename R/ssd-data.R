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

#' Data from fitdists Object
#'
#' Get a tibble of the original data.
#'
#' @inheritParams params
#'
#' @return A tibble of the original data.
#' @seealso [`augment.fitdists()`],  [`ssd_ecd_data()`] and [`ssd_sort_data()`]
#' @export
#'
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' ssd_data(fits)
ssd_data <- function(x) {
  chk_s3_class(x, "fitdists")
  .org_data_fitdists(x)
}
