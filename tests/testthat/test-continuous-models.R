cesa_continuous = load_cesa(get_test_file("cesa_hg38_for_test.rds"))


make_continuous_likelihood_inputs = function() {
  covariate_data = data.table::data.table(
    Unique_Patient_Identifier = c("S1", "S2", "S3", "S4"),
    covariate_value = c(20, 40, 60, 80)
  )

  list(
    covariate_data = covariate_data,
    rates_tumors_with = c(S1 = 0.010, S3 = 0.020),
    rates_tumors_without = c(S2 = 0.015, S4 = 0.005)
  )
}


manual_continuous_minusloglik = function(gamma_with,
                                         gamma_without,
                                         rates_tumors_with,
                                         rates_tumors_without) {
  z_with = gamma_with * rates_tumors_with
  z_without = gamma_without * rates_tumors_without

  -(sum(log(-expm1(-z_with))) - sum(z_without))
}


test_that("continuous likelihood functions use the expected parameterizations", {
  x = make_continuous_likelihood_inputs()

  linear_lik = sswm_age_lik(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  logistic_lik = sswm_age_lik_logistic(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  sigmoid_lik = sswm_age_lik_sigmoid(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  expect_identical(bbmle::parnames(linear_lik), c("beta0", "beta1"))
  expect_identical(bbmle::parnames(logistic_lik), c("L", "k", "m"))
  expect_identical(bbmle::parnames(sigmoid_lik), c("C", "L", "s", "m"))
})


test_that("linear continuous likelihood matches direct calculation", {
  x = make_continuous_likelihood_inputs()

  lik = sswm_age_lik(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  params = c(beta0 = 2, beta1 = 0.01)

  covariate_by_sample = setNames(
    x$covariate_data$covariate_value,
    x$covariate_data$Unique_Patient_Identifier
  )

  gamma_with = params["beta0"] +
    params["beta1"] * covariate_by_sample[names(x$rates_tumors_with)]

  gamma_without = params["beta0"] +
    params["beta1"] * covariate_by_sample[names(x$rates_tumors_without)]

  expected = manual_continuous_minusloglik(
    gamma_with = gamma_with,
    gamma_without = gamma_without,
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without
  )

  expect_equal(lik(params), expected, tolerance = 1e-12)
})


test_that("logistic continuous likelihood matches direct calculation", {
  x = make_continuous_likelihood_inputs()

  lik = sswm_age_lik_logistic(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  params = c(L = log(4), k = 0.08, m = 50)

  covariate_by_sample = setNames(
    x$covariate_data$covariate_value,
    x$covariate_data$Unique_Patient_Identifier
  )

  gamma = function(z) {
    exp(params["L"]) * plogis(params["k"] * (z - params["m"]))
  }

  expected = manual_continuous_minusloglik(
    gamma_with = gamma(covariate_by_sample[names(x$rates_tumors_with)]),
    gamma_without = gamma(covariate_by_sample[names(x$rates_tumors_without)]),
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without
  )

  expect_equal(lik(params), expected, tolerance = 1e-12)
})


test_that("generalized sigmoid likelihood matches direct calculation", {
  x = make_continuous_likelihood_inputs()

  lik = sswm_age_lik_sigmoid(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  params = c(C = log(0.8), L = log(4), s = 2, m = 50)

  covariate_by_sample = setNames(
    x$covariate_data$covariate_value,
    x$covariate_data$Unique_Patient_Identifier
  )

  gamma = function(z) {
    weight = plogis(params["s"] * (log(z) - log(params["m"])))
    exp(params["C"]) +
      (exp(params["L"]) - exp(params["C"])) * weight
  }

  expected = manual_continuous_minusloglik(
    gamma_with = gamma(covariate_by_sample[names(x$rates_tumors_with)]),
    gamma_without = gamma(covariate_by_sample[names(x$rates_tumors_without)]),
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without
  )

  expect_equal(lik(params), expected, tolerance = 1e-12)
})


test_that("continuous likelihoods match covariates by patient ID, not row order", {
  x = make_continuous_likelihood_inputs()
  shuffled_covariates = x$covariate_data[c(4, 1, 3, 2)]

  linear_a = sswm_age_lik(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )
  linear_b = sswm_age_lik(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = shuffled_covariates
  )

  logistic_a = sswm_age_lik_logistic(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )
  logistic_b = sswm_age_lik_logistic(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = shuffled_covariates
  )

  sigmoid_a = sswm_age_lik_sigmoid(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )
  sigmoid_b = sswm_age_lik_sigmoid(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = shuffled_covariates
  )

  expect_equal(
    linear_a(c(beta0 = 2, beta1 = 0.01)),
    linear_b(c(beta0 = 2, beta1 = 0.01)),
    tolerance = 1e-12
  )

  expect_equal(
    logistic_a(c(L = log(4), k = 0.08, m = 50)),
    logistic_b(c(L = log(4), k = 0.08, m = 50)),
    tolerance = 1e-12
  )

  expect_equal(
    sigmoid_a(c(C = log(0.8), L = log(4), s = 2, m = 50)),
    sigmoid_b(c(C = log(0.8), L = log(4), s = 2, m = 50)),
    tolerance = 1e-12
  )
})


test_that("invalid continuous-model parameters return a finite penalty", {
  x = make_continuous_likelihood_inputs()

  linear_lik = sswm_age_lik(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  sigmoid_lik = sswm_age_lik_sigmoid(
    rates_tumors_with = x$rates_tumors_with,
    rates_tumors_without = x$rates_tumors_without,
    covariate_data = x$covariate_data
  )

  bad_linear = linear_lik(c(beta0 = -1, beta1 = 0))
  bad_sigmoid = sigmoid_lik(c(C = 0, L = 1, s = 2, m = -1))

  expect_true(is.finite(bad_linear))
  expect_true(is.finite(bad_sigmoid))
  expect_gt(bad_linear, 1e50)
  expect_gt(bad_sigmoid, 1e50)
})


test_that("continuous ces_variant functions add model-specific results", {
  # FOXA1 F266L is already used in test-models.R and is recurrent in this fixture.
  cesa = copy_cesa(cesa_continuous)

  # Use a positive synthetic continuous covariate so it is valid for all
  # three models, including the generalized sigmoid model.
  cesa@samples[, continuous_test_covariate := seq_len(.N)]

  covariate_data = cesa@samples[
    ,
    .(
      Unique_Patient_Identifier,
      covariate_value = continuous_test_covariate
    )
  ]

  variant = select_variants(cesa, variant_ids = "FOXA1 F266L")
  expect_equal(nrow(variant), 1)

  invisible(capture.output({
    
    cesa = ces_variant_linear(
      cesa,
      variants = variant,
      model = sswm_age_lik,
      lik_args = list(covariate_data = covariate_data),
      optimizer = "COBYLA",
      hold_out_same_gene_samples = F,
      conf = NULL,
      run_name = "linear_test",
      cores = 1
    )
    
    cesa = ces_variant_logistic(
      cesa,
      variants = variant,
      model = sswm_age_lik_logistic,
      lik_args = list(covariate_data = covariate_data),
      optimizer = "COBYLA",
      hold_out_same_gene_samples = F,
      conf = NULL,
      run_name = "logistic_test",
      cores = 1
    )
    
    cesa = ces_variant_sigmoid(
      cesa,
      variants = variant,
      model = sswm_age_lik_sigmoid,
      lik_args = list(covariate_data = covariate_data),
      optimizer = "COBYLA",
      hold_out_same_gene_samples = F,
      conf = NULL,
      run_name = "sigmoid_test",
      cores = 1
    )
    
  }))

  expect_true(all(
    c("linear_test", "logistic_test", "sigmoid_test") %in%
      names(cesa@selection_results)
  ))

  linear_result = cesa@selection_results$linear_test
  logistic_result = cesa@selection_results$logistic_test
  sigmoid_result = cesa@selection_results$sigmoid_test

  expect_equal(nrow(linear_result), 1)
  expect_equal(nrow(logistic_result), 1)
  expect_equal(nrow(sigmoid_result), 1)

  expect_true(all(
    c("variant_name", "variant_id", "beta0", "beta1", "loglikelihood") %in%
      names(linear_result)
  ))

  expect_true(all(
    c("variant_name", "variant_id", "L", "k", "m", "loglikelihood") %in%
      names(logistic_result)
  ))

  expect_true(all(
    c("variant_name", "variant_id", "C", "L", "s", "m", "loglikelihood") %in%
      names(sigmoid_result)
  ))

  expect_true(all(is.finite(
    unlist(linear_result[, .(beta0, beta1, loglikelihood)])
  )))

  expect_true(all(is.finite(
    unlist(logistic_result[, .(L, k, m, loglikelihood)])
  )))

  expect_true(all(is.finite(
    unlist(sigmoid_result[, .(C, L, s, m, loglikelihood)])
  )))

  # Linear selection intensity is constrained to be positive over the
  # observed covariate range.
  x_range = range(covariate_data$covariate_value)
  gamma_endpoints =
    linear_result$beta0 + linear_result$beta1 * x_range
  expect_true(all(gamma_endpoints >= 1e-3 - 1e-6))

  # Smoke-test the plotting function using the fitted results.
  p = plot_effects_continuous(
    cesa,
    linear_run_name = "linear_test",
    logistic_run_name = "logistic_test",
    sigmoid_run_name = "sigmoid_test",
    variants = linear_result$variant_name,
    covariate_col = "continuous_test_covariate",
    output = "continuous",
    continuous_models = "all",
    x_title = "Test covariate"
  )

  expect_true(inherits(p, "ggplot"))
})
