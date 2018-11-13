cmake_minimum_required (VERSION 3.10)
project (FESR)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Boost)
find_package(GSL REQUIRED)
find_package(ROOT COMPONENTS MINUIT)

# includes ROOT header files
include(${ROOT_USE_FILE})

set(DIR src)
include_directories(${DIR} ${BOOST_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
# file(GLOB SRC "*.cpp")
file(GLOB SOURCES src/*.cpp)

add_executable(FESR main.cpp ${SOURCES})
target_link_libraries(FESR GSL::gsl ${ROOT_LIBRARIES})

#################################
# Testing
#################################
if(test)
  # ctest
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

  set (TEST_DIR test)
  file(GLOB TEST_SOURCES test/*cpp)

  #################
  # Unit Tests
  #################
  add_executable(runUnitTests ${SOURCES} ${TEST_SOURCES})

  # Standard linking to gtest stuff.
  target_link_libraries(runUnitTests GTest::GTest GTest::Main GSL::gsl)

  # Extra linking for the project.
  # target_link_libraries(...)

  add_test(AllTestsInFESR runUnitTests)
endif()
