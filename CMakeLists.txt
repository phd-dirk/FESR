cmake_minimum_required (VERSION 3.10)
project (FESR)
set (CMAKE_CXX_STANDARD 17)

find_package(Boost 1.66.0)

set(DIR src)
include_directories(${DIR} ${BOOST_INCLUDE_DIRS})
set(SOURCES main.cpp ${DIR}/data.hpp)

add_executable(FESR ${SOURCES})
target_link_libraries(FESR)

#################################
# Testing
#################################
if(test)
  # ctest
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS})

  set (TEST_DIR test)
  set (TEST_SOURCES ${TEST_DIR}/constants_test.cpp ${TEST_DIR}/weights_test.cpp ${TEST_DIR}/data_test.cpp)

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
