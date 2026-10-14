cmake_minimum_required(VERSION 3.16...3.27)
project(xenonCryostatSim LANGUAGES CXX)

# Requires a Geant4 (>= 11.0) built with:
#   -DGEANT4_USE_GDML=ON        (Xerces-C; needed for the alpha-source import)
#   -DGEANT4_USE_QT=ON / OpenGL (only if you want the interactive viewer)
#   -DGEANT4_BUILD_MULTITHREADED=ON (recommended; analysis merging is enabled)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# Placeholder geometry guard.
#   ON  (default): GeometryParameters.hh GEO_V(...) placeholders are used; the
#                  run prints a DEMO banner and writes xenon_sim_DEMO.root.
#   OFF          : unreplaced placeholders are NaN and the app refuses to build
#                  the geometry, listing what is missing.  Use for CAD runs.
option(GEO_DEMO "Use placeholder DEMO geometry values (NOT the CAD numbers)" ON)
if(GEO_DEMO)
  add_compile_definitions(GEO_USE_DEMO_VALUES)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/generated)

# Git version for the run provenance: generated at BUILD time (not just at configure time), so the
# commit recorded in the output is the one that was compiled. See cmake/GitVersion.cmake.
set(GIT_VERSION_HEADER ${PROJECT_BINARY_DIR}/generated/GitVersion.hh)
add_custom_target(git_version ALL
  COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR} -DOUT=${GIT_VERSION_HEADER}
          -P ${PROJECT_SOURCE_DIR}/cmake/GitVersion.cmake
  BYPRODUCTS ${GIT_VERSION_HEADER})

file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

# geometry_probe.cc has its own main(): keep it out of the application.
list(FILTER sources EXCLUDE REGEX "geometry_probe\\.cc$")

add_executable(xenonSim xenonSim.cc ${sources} ${headers})
target_link_libraries(xenonSim ${Geant4_LIBRARIES})
add_dependencies(xenonSim git_version)

# Geometry self-check: validates the header and probes the solids (no physics,
# no data files needed).  Run it every time you change GeometryParameters.hh.
add_executable(geometryProbe src/geometry_probe.cc ${headers})
target_link_libraries(geometryProbe ${Geant4_LIBRARIES})

# Stage runtime assets so ./xenonSim macros/run.mac works from the build dir.
file(GLOB macro_files CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/macros/*.mac)
foreach(_a ${macro_files})
  get_filename_component(_n ${_a} NAME)
  configure_file(${_a} ${PROJECT_BINARY_DIR}/macros/${_n} COPYONLY)
endforeach()
configure_file(${PROJECT_SOURCE_DIR}/gdml/alpha_source.gdml
               ${PROJECT_BINARY_DIR}/gdml/alpha_source.gdml COPYONLY)
