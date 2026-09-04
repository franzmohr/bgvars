test_that("create_gvecmodel returns the documented structure", {
  submodel_data <- gvec_data()
  object <- create_gvecmodel(submodel_data = submodel_data,
                             global_data = gvec_global_data())

  expect_s3_class(object, "gvecmodel")
  expect_false(inherits(object, "gvarmodel"))
  expect_equal(names(object), c("global", "weights", "submodels"))
  expect_null(object[["weights"]])
  expect_null(object[["submodels"]])
})

test_that("create_gvecmodel builds the same global index as create_gvarmodel", {
  submodel_data <- gvec_data()

  gvec <- create_gvecmodel(submodel_data = submodel_data,
                           global_data = gvec_global_data())
  gvar <- create_gvarmodel(submodel_data = submodel_data,
                           global_data = gvec_global_data())

  expect_equal(gvec[["global"]][["index"]], gvar[["global"]][["index"]])
  expect_equal(gvec[["global"]][["endogen"]], gvar[["global"]][["endogen"]])
})

test_that("create_gvecmodel validates its arguments", {
  submodel_data <- gvec_data()

  expect_error(create_gvecmodel(submodel_data = unclass(submodel_data)),
               "must be of class 'submodeldata'")
  expect_error(create_gvecmodel(submodel_data = submodel_data,
                                global_data = as.data.frame(gvec_global_data())),
               "must be of class 'ts'")
})
