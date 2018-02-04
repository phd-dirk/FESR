cmake_minimum_required (VERSION 3.10)
project (FESR)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Boost 1.66.0)
find_package(GSL REQUIRED)

set(DIR src lib/cRunDec)
include_directories(${DIR} ${BOOST_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
set(SOURCES main.cpp lib/cRunDec/CRunDec.cpp)

add_executable(FESR ${SOURCES})
target_link_libraries(FESR GSL::gsl)

#################################
# Testing
#################################
if(test)
  # ctest
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS})

  set (TEST_DIR test)
  set (TEST_SOURCES ${TEST_DIR}/constants_test.cpp ${TEST_DIR}/weights_test.cpp
    ${TEST_DIR}/data_test.cpp)

  #################
  # Unit Tests
  #################
  add_executable(runUnitTests ${TEST_SOURCES})

  # Standard linking to gtest stuff.
  target_link_libraries(runUnitTests GTest::GTest GTest::Main)

  # Extra linking for the project.
  # target_link_libraries(...)

  add_test(AllTestsInFESR runUnitTests)
endif()
