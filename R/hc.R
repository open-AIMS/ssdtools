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

#' Hazard Concentrations for Species Sensitivity Distributions
#'
#' Gets concentration(s) that protect specified percentage(s) of species.
#'
#' If `ci = TRUE` and `parametric = TRUE` uses parameteric bootstrapping to get confidence intervals on the
#' hazard concentrations(s).  If `ci = TRUE` and `parametric = FALSE` uses non-parameteric bootstrapping to get confidence intervals on the
#' hazard concentrations(s).
#'
#' @inheritParams params
#' @return A tibble of corresponding hazard concentrations.
#' @seealso [`predict.fitdists()`] and [`ssd_hp()`].
#' @export
#' @examples
#' fits <- ssd_fit_dists(ssddata::ccme_boron)
#' ssd_hc(fits)
#' ssd_hc(estimates(fits))
#' ssd_hc(ssd_match_moments())
ssd_hc <- function(x, ...) {
  UseMethod("ssd_hc")
}

no_ssd_hc <- function() {
  tibble(
    dist = character(0),
    percent = numeric(0),
    est = numeric(0),
    se = numeric(0),
    lcl = numeric(0),
    ucl = numeric(0),
    wt = numeric(0),
    nboot = integer(0),
    pboot = numeric(0)
  )
}

.ssd_hc_dist <- function(x, dist, proportion) {
  fun <- paste0("ssd_q", dist)
  args <- list(p = proportion)
  args <- c(as.list(x), args)
  est <- do.call(fun, args)
  tibble(
    dist = dist,
    percent = proportion * 100, est = est,
    se = NA_real_, lcl = NA_real_, ucl = NA_real_,
    wt = 1,
    nboot = 0L, pboot = NA_real_
  )
}

.ssd_hc_tmbfit <- function(x, proportion, ci, level, nboot, min_pboot,
                           data, rescale, weighted, censoring, min_pmix,
                           range_shape1, range_shape2, parametric, control) {
  args <- estimates(x)
  args$p <- proportion
  dist <- .dist_tmbfit(x)
  what <- paste0("ssd_q", dist)
  est <- do.call(what, args)
  if (!ci) {
    na <- rep(NA_real_, length(proportion))
    return(tibble(
      dist = rep(dist, length(proportion)),
      percent = proportion * 100,
      est = est * rescale,
      se = na,
      lcl = na,
      ucl = na,
      wt = rep(1, length(proportion)),
      nboot = rep(0L, length(proportion)),
      pboot = na,
      samples = na
    ))
  }
  censoring <- censoring / rescale
  fun <- safely(fit_tmb)
  estimates <- boot_estimates(x,
    fun = fun, nboot = nboot, data = data, weighted = weighted,
    censoring = censoring, min_pmix = min_pmix,
    range_shape1 = range_shape1,
    range_shape2 = range_shape2,
    parametric = parametric,
    control = control
  )
  
  samples <- sample_estimates(estimates, what, x = proportion)
  cis <- cis_estimates(estimates, what, level = level, x = proportion)
  hc <- tibble(
    dist = dist,
    percent = proportion * 100, est = est * rescale,
    se = cis$se * rescale, lcl = cis$lcl * rescale, ucl = cis$ucl * rescale,
    wt = rep(1, length(proportion)),
    nboot = nboot, pboot = length(estimates) / nboot,
    samples = samples
  )
  replace_min_pboot_na(hc, min_pboot)
}

.ssd_hc_fitdists <- function(x, percent, ci, level, nboot,
                             average, min_pboot, parametric, control) {
  if (!length(x) || !length(percent)) {
    return(no_ssd_hc())
  }

  if (is.null(control)) {
    control <- .control_fitdists(x)
  }
  data <- .data_fitdists(x)
  rescale <- .rescale_fitdists(x)
  censoring <- .censoring_fitdists(x)
  min_pmix <- .min_pmix_fitdists(x)
  range_shape1 <- .range_shape1_fitdists(x)
  range_shape2 <- .range_shape2_fitdists(x)
  weighted <- .weighted_fitdists(x)
  unequal <- .unequal_fitdists(x)
  weight <- glance(x)$weight
  names(weight) <- glance(x)$dist
  
  if (parametric && ci && identical(censoring, c(NA_real_, NA_real_))) {
    wrn("Parametric CIs cannot be calculated for inconsistently censored data.")
    ci <- FALSE
  }

  if (parametric && ci && unequal) {
    wrn("Parametric CIs cannot be calculated for unequally weighted data.")
    ci <- FALSE
  }
  if (!ci) {
    nboot <- 0L
  }
  seeds <- seed_streams(length(x))

  if (!average) {
    hc <- future_map(x, .ssd_hc_tmbfit,  
      proportion = percent / 100, ci = ci, level = level, nboot = nboot,
      min_pboot = min_pboot,
      data = data, rescale = rescale, weighted = weighted, 
      censoring = censoring,
      min_pmix = min_pmix, range_shape1 = range_shape1, 
      range_shape2 = range_shape2,
      parametric = parametric, control = control,
      .options = furrr::furrr_options(seed = seeds)
    )    
    hc <- mapply(
      function(x, y) {
        x$wt <- y
        x
      },
      x = hc, y = weight,
      USE.NAMES = FALSE, SIMPLIFY = FALSE
    )
    hc <- bind_rows(hc)
    hc$method <- if (parametric) "parametric" else "non-parametric"
    hc <- hc[c("dist", "percent", "est", "se", "lcl", "ucl", "wt", 
               "method", "nboot", "pboot")]
    return(hc)
  }
 
  nboot_vals <- round(round(nboot*weight))
  # # check values are 1 or greater
  # nboot_valid <- names(nboot_vals)[which(nboot_vals>0)]
  # 
  # # remove distributions where nboot==0
  # x <- x[nboot_valid]
  # nboot_vals[nboot_valid]
  
  hc <- furrr::future_map2(.x = x, .y = nboot_vals,  
       ~ .ssd_hc_tmbfit(x = .x, proportion = percent / 100, ci = ci, 
                        level = level, 
                        nboot = .y,
                        min_pboot = min_pboot,
                        data = data, rescale = rescale, weighted = weighted, 
                        censoring = censoring,
                        min_pmix = min_pmix, range_shape1 = range_shape1, 
                        range_shape2 = range_shape2,
                        parametric = parametric, control = control),  
       .options = furrr::furrr_options(seed = seeds)
  ) 

  hc <- lapply(hc, FUN = function(x) x |> tidyr::unnest_longer(samples)) |> 
    bind_rows()
  pboot_chk <- hc |> dplyr::select(dist, pboot) |> unique() |> 
    dplyr::filter(pboot<min_pboot)
  dists_fail <- paste(pboot_chk$dist, collapse = "; ")

  if(nrow(pboot_chk)>0) {
    stop(paste("The ", dists_fail, " distribution(s) fail(s) the minimum bootstrap convergence criteria of ", 
               min_pboot, ". Please drop the failing distribution(s), or modify pboot.", sep=""))
  }
  
  new_pboot <- nrow(hc)/length(percent)/nboot  
  method <- if (parametric) "parametric" else "non-parametric"    
  
  hc |> dplyr::select(percent, samples) |> 
    dplyr::group_by(percent) |> 
    dplyr::summarise(est = mean(samples),
                     lcl = quantile(samples, probs = probs(level)[1]),
                     ucl = quantile(samples, probs = probs(level)[2]),
                     se = sd(samples)) |> 
    dplyr::mutate(dist = "average",
                  method = method,
                  nboot = nboot, 
                  pboot = new_pboot,
                  wt = 1) |> 
    dplyr::select(dist, percent, est, se, lcl, ucl, wt, method, nboot, pboot)
  
}

#' @describeIn ssd_hc Hazard Concentrations for Distributional Estimates
#' @export
ssd_hc.list <- function(x, percent = 5, ...) {
  chk_list(x)
  chk_named(x)
  chk_unique(names(x))
  chk_unused(...)

  if (!length(x)) {
    return(no_ssd_hc())
  }
  hc <- mapply(.ssd_hc_dist, x, names(x),
    MoreArgs = list(proportion = percent / 100),
    SIMPLIFY = FALSE
  )
  bind_rows(hc)
}

#' @describeIn ssd_hc Hazard Concentrations for fitdists Object
#' @export
ssd_hc.fitdists <- function(x, percent = 5, ci = FALSE, level = 0.95, nboot = 1000,
                            average = TRUE, delta = 7, min_pboot = 0.99,
                            parametric = TRUE,
                            control = NULL, ...) {
  chk_vector(percent)
  chk_numeric(percent)
  chk_range(percent, c(0, 100))
  chk_flag(ci)
  chk_number(level)
  chk_range(level)
  chk_whole_number(nboot)
  chk_gt(nboot)
  chk_flag(average)
  chk_number(delta)
  chk_gte(delta)
  chk_number(min_pboot)
  chk_range(min_pboot)
  chk_flag(parametric)
  chk_null_or(control, vld = vld_list)
  chk_unused(...)

  x <- subset(x, delta = delta)
  hc <- .ssd_hc_fitdists(x, percent,
    ci = ci, level = level, nboot = nboot, min_pboot = min_pboot, control = control,
    average = average, parametric = parametric
  )
  warn_min_pboot(hc, min_pboot)
}

#' @describeIn ssd_hc Hazard Concentrations for fitburrlioz Object
#' '
#' @export
#' @examples
#' fit <- ssd_fit_burrlioz(ssddata::ccme_boron)
#' ssd_hc(fit)
#'
#' @export
ssd_hc.fitburrlioz <- function(x, percent = 5, ci = FALSE, level = 0.95, nboot = 1000,
                               min_pboot = 0.99, parametric = FALSE, ...) {
  check_dim(x, values = 1L)
  chk_named(x)
  chk_subset(names(x), c("burrIII3", "invpareto", "llogis", "lgumbel"))
  chk_vector(percent)
  chk_numeric(percent)
  chk_range(percent, c(0, 100))
  chk_flag(ci)
  chk_number(level)
  chk_range(level)
  chk_whole_number(nboot)
  chk_gt(nboot)
  chk_number(min_pboot)
  chk_range(min_pboot)
  chk_flag(parametric)
  chk_unused(...)

  if (names(x) != "burrIII3" || !ci || !length(percent)) {
    class(x) <- class(x)[-1]
    return(ssd_hc(x,
      percent = percent, ci = ci, level = level,
      nboot = nboot, min_pboot = min_pboot,
      average = FALSE, parametric = parametric
    ))
  }
  hc <- .ssd_hc_burrlioz_fitdists(x,
    percent = percent, level = level, nboot = nboot,
    min_pboot = min_pboot, parametric = parametric
  )
  warn_min_pboot(hc, min_pboot)
}
