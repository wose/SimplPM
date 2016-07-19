set(test_SRC
  main.cpp
  )

add_executable(test_suite ${test_SRC})
target_link_libraries(test_suite
  gmock
  gtest
  simplpm
  )

add_test(NAME Test
  COMMAND test_suite)
