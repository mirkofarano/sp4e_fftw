cmake_minimum_required (VERSION 3.1)
project (Particles)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

################################################################
# libpart
################################################################

add_library(part
  compute_boundary.cc
  compute_verlet_integration.cc 
  particle.cc 
  planet.cc
  compute_gravity.cc 
  csv_reader.cc 
  particles_factory_interface.cc 
  planets_factory.cc 
  compute_contact.cc 
  compute_kinetic_energy.cc 
  csv_writer.cc 
  system.cc 
  compute_energy.cc 
  compute_potential_energy.cc 
  ping_pong_ball.cc
  material_point.cc 
  system_evolution.cc 
  ping_pong_balls_factory.cc 
  compute_interaction.cc
  compute_temperature.cc
  material_points_factory.cc 
  )

add_executable(particles main.cc)
target_link_libraries(particles part)

################################################################
# Google test
################################################################


include(GoogleTest)
enable_testing()
find_package(GTest)
if (GTEST_FOUND)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_kepler test_kepler.cc)
add_executable(test_fft test_fft.cc)
target_link_libraries(test_kepler part ${GTEST_BOTH_LIBRARIES} pthread)
target_link_libraries(test_fft part ${GTEST_BOTH_LIBRARIES} ${FFTW_LIBRARIES} pthread)
gtest_discover_tests(test_kepler)
gtest_discover_tests(test_fft)
endif()

################################################################
# Doxygen
################################################################

find_package(Doxygen)
if (DOXYGEN_FOUND)

# to set other options, read: https://cmake.org/cmake/help/v3.9/module/FindDoxygen.html
doxygen_add_docs(
    doxygen
    ${PROJECT_SOURCE_DIR}
    COMMENT "Generate html pages"
)
add_custom_target(doc DEPENDS doxygen)
endif(DOXYGEN_FOUND)
