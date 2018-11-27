test_that("Test that optimal design with point prior runs through", {

  opt_design(alpha = .025,
             beta = .2,
             beta.2 = .2,
             lambda = c(0, 0, 0),
             prior = "point",
             power = "equal",
             delta.mcr = .1,
             delta.alt = .3,
             tau = .1,
             n.max = Inf,
             distribution = "t",
             n1 = NA)

  expect_equal(1,1)

})



test_that("Test that optimal design with normal prior runs through", {

  opt_design(alpha = .025,
             beta = .2,
             beta.2 = .2,
             lambda = c(1, 1, 1),
             prior = "normal",
             power = "predictive",
             delta.mcr = .1,
             delta.alt = .3,
             tau = .1,
             n.max = 150,
             distribution = "normal",
             n1 = 60)

  expect_equal(1,1)

})

