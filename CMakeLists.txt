cmake_minimum_required (VERSION 3.10)
project (FESR)
set (CMAKE_CXX_STANDARD 14)

add_executable(FESR main.cpp)
target_link_libraries(FESR)

#################################
# Testing
#################################
if(test)
  # ctest
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS})

  #################
  # Unit Tests
  #################
  add_executable(runUnitTests test.cpp)

  # Standard linking to gtest stuff.
  target_link_libraries(runUnitTests GTest::GTest GTest::Main)

  # Extra linking for the project.
  # target_link_libraries(...)

  add_test(AllTestsInFESR runUnitTests)
endif()
