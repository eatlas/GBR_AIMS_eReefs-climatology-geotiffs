source("time-functions.R", local = TRUE)

# Test suite for convert_netcdf_time function
test_netcdf_time_conversion <- function() {
  # Create a simple test helper function
  run_test <- function(desc, time_values, time_units, expected_results) {
    cat(sprintf("\nTest: %s\n", desc))
    
    # Use tryCatch to handle potential errors
    result <- tryCatch({
      convert_netcdf_time(time_values, time_units)
    }, error = function(e) {
      cat(sprintf("  Error: %s\n", e$message))
      return(NA)
    })
    
    # Check for NAs in result
    if (any(is.na(result))) {
      cat("  Result: FAILED (NA values in result)\n")
      return(FALSE)
    }
    
    passed <- all(format(result) == format(expected_results))
    cat(sprintf("  Result: %s\n", ifelse(passed, "PASSED", "FAILED")))
    if (!passed) {
      cat("  Expected: ", format(expected_results), "\n")
      cat("  Got:      ", format(result), "\n")
    }
    return(passed)
  }

  # Function to test for expected errors
  run_error_test <- function(desc, time_values, time_units, expected_error_pattern) {
    cat(sprintf("\nTest: %s\n", desc))
    
    # Use tryCatch to verify an error is thrown
    result <- tryCatch({
      convert_netcdf_time(time_values, time_units)
      cat("  Result: FAILED (No error was thrown, but one was expected)\n")
      return(FALSE)
    }, error = function(e) {
      # Check if error message matches expected pattern
      if (grepl(expected_error_pattern, e$message)) {
        cat(sprintf("  Result: PASSED (Got expected error: %s)\n", 
                    substr(e$message, 1, 60)))
        return(TRUE)
      } else {
        cat(sprintf("  Result: FAILED (Error message didn't match pattern)\n"))
        cat(sprintf("  Expected pattern: %s\n", expected_error_pattern))
        cat(sprintf("  Got error: %s\n", e$message))
        return(FALSE)
      }
    })
    
    return(result)
  }
  
  # Test 1: Basic days format with +10 timezone
  test1 <- run_test(
    "Basic days format with +10 timezone",
    c(0, 365, 730),
    "days since 1990-01-01 00:00:00 +10",
    c(as.POSIXct("1989-12-31 14:00:00", tz="UTC"),
      as.POSIXct("1990-12-31 14:00:00", tz="UTC"),
      as.POSIXct("1991-12-31 14:00:00", tz="UTC"))
  )
  
  # Test 2: Days format with no timezone (should assume UTC)
  test2 <- run_test(
    "Days format with no timezone",
    c(0, 365, 730),
    "days since 1990-01-01 00:00:00",
    c(as.POSIXct("1990-01-01 00:00:00", tz="UTC"),
      as.POSIXct("1991-01-01 00:00:00", tz="UTC"),
      as.POSIXct("1992-01-01 00:00:00", tz="UTC"))
  )
  
  # Test 3: Seconds format
  test3 <- run_test(
    "Seconds format",
    c(0, 86400, 172800),
    "seconds since 1990-01-01 00:00:00",
    c(as.POSIXct("1990-01-01 00:00:00", tz="UTC"),
      as.POSIXct("1990-01-02 00:00:00", tz="UTC"),
      as.POSIXct("1990-01-03 00:00:00", tz="UTC"))
  )
  
  # Test 4: Hours format with negative timezone
  test4 <- run_test(
    "Hours format with negative timezone",
    c(0, 24, 48),
    "hours since 1990-01-01 00:00:00 -05",
    c(as.POSIXct("1990-01-01 05:00:00", tz="UTC"),
      as.POSIXct("1990-01-02 05:00:00", tz="UTC"),
      as.POSIXct("1990-01-03 05:00:00", tz="UTC"))
  )
  
  # Test 5: Full format timezone (+1000)
  test5 <- run_test(
    "Full format timezone (+1000)",
    c(0, 365),
    "days since 1990-01-01 00:00:00 +1000",
    c(as.POSIXct("1989-12-31 14:00:00", tz="UTC"),
      as.POSIXct("1990-12-31 14:00:00", tz="UTC"))
  )
  
  # Test 6: Leap year handling
  test6 <- run_test(
    "Leap year handling",
    c(59, 60, 61),  # Feb 28, 29, Mar 1 in 2020
    "days since 2020-01-01 00:00:00",
    c(as.POSIXct("2020-02-29 00:00:00", tz="UTC"),
      as.POSIXct("2020-03-01 00:00:00", tz="UTC"),
      as.POSIXct("2020-03-02 00:00:00", tz="UTC"))
  )
  
  # Test 7: Minutes format
  test7 <- run_test(
    "Minutes format",
    c(0, 60, 120),
    "minutes since 1990-01-01 00:00:00",
    c(as.POSIXct("1990-01-01 00:00:00", tz="UTC"),
      as.POSIXct("1990-01-01 01:00:00", tz="UTC"),
      as.POSIXct("1990-01-01 02:00:00", tz="UTC"))
  )
  
  # Test 8: Edge case with very large time values
  test8 <- run_test(
    "Edge case with very large time values",
    c(36524),  # ~100 years
    "days since 1900-01-01 00:00:00",
    c(as.POSIXct("2000-01-01 00:00:00", tz="UTC"))
  )
  
  # Test 9: Single-digit timezone (+1)
  test9 <- run_test(
    "Single-digit timezone (+1)",
    c(0),
    "days since 1990-01-01 00:00:00 +1",
    c(as.POSIXct("1989-12-31 23:00:00", tz="UTC"))
  )
  # Test 10: Unsupported time unit
  test10 <- run_error_test(
    "Unsupported time unit",
    c(0, 1, 2),
    "months since 1990-01-01 00:00:00",
    "Unsupported time unit"  # Error message should contain this pattern
  )
  
  # Test 11: Malformed time unit string
  test11 <- run_error_test(
    "Malformed time unit string",
    c(0, 1, 2),
    "this is not a valid time unit string",
    "Unsupported time unit"  # Error message should contain this pattern
  )
  
  # Update overall results check to include new tests
  all_passed <- all(c(test1, test2, test3, test4, test5, test6, test7, test8, test9, 
                      test10, test11))
  cat(sprintf("\nOverall test result: %s\n", ifelse(all_passed, "ALL TESTS PASSED", "SOME TESTS FAILED")))
}

# Run the tests
test_netcdf_time_conversion()