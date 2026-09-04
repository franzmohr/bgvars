# Minimal stand-in for the output of bvartools::selection_criteria, reduced to
# the 'summary' element that choose_best_model works on.
selcrit_fixture <- function() {
  criteria <- list(
    list(summary = data.frame(LL = -10, AIC = 30, BIC = 40, HQ = 35)),
    list(summary = data.frame(LL = -8,  AIC = 20, BIC = 50, HQ = 33)),
    list(summary = data.frame(LL = -9,  AIC = 25, BIC = 45, HQ = 31)),
    list(summary = data.frame(LL = -20, AIC = 60, BIC = 70, HQ = 65)),
    list(summary = data.frame(LL = -30, AIC = 55, BIC = 65, HQ = 60))
  )
  names(criteria) <- c("US", "US", "US", "JP", "JP")
  class(criteria) <- c("submodelselcritlist", "list")
  criteria
}

test_that("choose_best_model minimises the information criterion", {
  criteria <- selcrit_fixture()

  # BIC is smallest for the first US model and the second JP model.
  expect_equal(choose_best_model(criteria), c(1, 5))
  expect_equal(choose_best_model(criteria, criterion = "BIC"), c(1, 5))
  expect_equal(choose_best_model(criteria, criterion = "AIC"), c(2, 5))
  expect_equal(choose_best_model(criteria, criterion = "HQ"), c(3, 5))
})

test_that("choose_best_model maximises the log-likelihood", {
  criteria <- selcrit_fixture()

  expect_equal(choose_best_model(criteria, criterion = "LL"), c(2, 4))
})

test_that("choose_best_model returns one model per sub-model", {
  criteria <- selcrit_fixture()
  best <- choose_best_model(criteria)

  expect_setequal(names(criteria)[best], unique(names(criteria)))
  expect_length(best, length(unique(names(criteria))))
})

# Minimal stand-in for a list of estimated sub-models.
submodelest_fixture <- function() {
  object <- list("a", "b", "c", "d")
  names(object) <- c("US", "US", "JP", "JP")
  class(object) <- c("submodelestlist", "list")
  object
}

test_that("select_submodels keeps the selected models and their names", {
  object <- submodelest_fixture()
  result <- select_submodels(object, c(2, 3))

  expect_s3_class(result, "submodelestlist")
  expect_equal(names(result), c("US", "JP"))
  expect_equal(result[[1]], "b")
  expect_equal(result[[2]], "c")
})

test_that("select_submodels works together with choose_best_model", {
  criteria <- selcrit_fixture()
  object <- list("m1", "m2", "m3", "m4", "m5")
  names(object) <- names(criteria)
  class(object) <- c("submodelestlist", "list")

  result <- select_submodels(object, choose_best_model(criteria))
  expect_equal(names(result), c("US", "JP"))
  expect_equal(unlist(result, use.names = FALSE), c("m1", "m5"))
})

test_that("select_submodels rejects selections that break the global model", {
  object <- submodelest_fixture()

  expect_error(select_submodels(unclass(object), 1:2),
               "must be of class 'submodelestlist'")
  expect_error(select_submodels(object, c(1, 2)),
               "more than one entity or country")
})
