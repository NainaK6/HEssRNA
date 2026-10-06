test_that("prcesDF4modelInhouse processes the data frame correctly", {

  # Example input data frame
  df4modelInhouseFmt <- data.frame(
    ID = c("A", "B"),
    R1 = c(1.23456, 2.34567),
    R2 = c(3.45678, 4.56789)
  )

  # Expected output after processing
  expected_output <- data.frame(
    ID = rep(c("A", "B"), each = 2),
    NoOfReplicates = c(1, 2, 1, 2),
    pwr = round(c(1.23456, 3.45678, 2.34567, 4.56789), 3) + 0.001
  )

  # Sort the expected output by ID and NoOfReplicates for consistent comparison
  expected_output <- expected_output[order(expected_output$ID, expected_output$NoOfReplicates), ]

  # Run the function
  output <- prcesDF4modelInhouse(df4modelInhouseFmt)

  # Sort the output similarly by ID and NoOfReplicates
  output <- output[order(output$ID, output$NoOfReplicates), ]

  # Test that the output has the expected structure and values
  expect_equal(output, expected_output)

})
