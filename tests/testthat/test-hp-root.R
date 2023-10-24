#    Copyright 2023 Australian Government Department of 
#    Climate Change, Energy, the Environment and Water
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

#    Copyright 2021 Province of British Columbia
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

test_that("hp root lnorm", {
  skip_on_os("linux") # FIXME
  fits <- ssd_fit_dists(ssddata::ccme_boron, dists = "lnorm")
  set.seed(102)
  hp_dist <- ssd_hp(fits, average = FALSE)
  hp_average <- ssd_hp(fits, average = TRUE)
  hp_root <- ssd_hp(fits, average = TRUE, root = TRUE)
  expect_identical(hp_average$est, hp_dist$est)
  expect_equal(hp_root, hp_average, tolerance = 1e-3)
  expect_equal(hp_average$est, 1.9543030195088, tolerance = 1e-6)
  expect_equal(hp_root$est, 1.95476926846743, tolerance = 1e-6)
  
  testthat::expect_snapshot({
    hp_root
  })
})

test_that("hp root all", {
  skip_on_os("linux") 
  fits <- ssd_fit_dists(ssddata::ccme_boron)
  set.seed(102)
  hp_average <- ssd_hp(fits, average = TRUE)
  hp_root <- ssd_hp(fits, average = TRUE, root = TRUE)
  expect_equal(hp_root, hp_average, tolerance = 1e-2)
  expect_equal(hp_average$est, 3.89879358571718, tolerance = 1e-6)
  expect_equal(hp_root$est, 3.91103597328257, tolerance = 1e-6)
  testthat::expect_snapshot({
    hp_root
  })
})

test_that("hp is hc", {
  skip_on_os("linux") 
  fits <- ssd_fit_dists(ssddata::ccme_boron)
  set.seed(102)
  conc <- 1
  hp_root <- ssd_hp(fits, conc = 1, average = TRUE, root = TRUE)
  hc_root <- ssd_hc(fits, percent = hp_root$est, average = TRUE, root = TRUE)
  expect_equal(hc_root$est, conc, tolerance = 1e-2)
  for(i in 1:100) {
    hp_root <- ssd_hp(fits, conc = hc_root$est, average = TRUE, root = TRUE)
    hc_root <- ssd_hc(fits, percent = hp_root$est, average = TRUE, root = TRUE)
  }
  skip("uniroot is biased...")
  expect_equal(hc_root$est, conc, tolerance = 1e-2)
})