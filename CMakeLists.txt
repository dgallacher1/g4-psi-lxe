cmake_minimum_required(VERSION 3.16...3.27)
project(xenonCryostatSim LANGUAGES CXX)

# Requires a Geant4 (>= 11.0) built with:
#   -DGEANT4_USE_GDML=ON        (Xerces-C; needed for the CAD/GDML alpha source)
#   -DGEANT4_USE_QT=ON / OpenGL (only if you want the interactive viewer)
#   -DGEANT4_BUILD_MULTITHREADED=ON (recommended; analysis merging is enabled below)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(xenonSim xenonSim.cc ${sources} ${headers})
target_link_libraries(xenonSim ${Geant4_LIBRARIES})

# Stage runtime assets next to the binary so ./xenonSim finds them.
set(ASSETS
  macros/beam.mac
  gdml/alpha_source.gdml
)
foreach(_a ${ASSETS})
  configure_file(${PROJECT_SOURCE_DIR}/${_a} ${PROJECT_BINARY_DIR}/${_a} COPYONLY)
endforeach()
