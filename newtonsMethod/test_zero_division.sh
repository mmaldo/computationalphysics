#!/bin/bash

# Compile the test runner with the isolated function
gcc newtonsMethod/fun_for_testing.c newtonsMethod/test_runner.c -o newtonsMethod/test_executable

# Run the executable and capture stderr
output=$(./newtonsMethod/test_executable 2>&1)
exit_code=$?

# Clean up the executable
rm newtonsMethod/test_executable

# Check for the non-zero exit code and the specific error message
if [ $exit_code -ne 0 ] && [[ "$output" == *"Error: Division by zero in fun()."* ]]; then
  echo "Test passed: Program exited with a non-zero status code and the correct error message."
  exit 0
else
  echo "Test failed: Program did not exit with a non-zero status code or the correct error message was not printed."
  exit 1
fi
