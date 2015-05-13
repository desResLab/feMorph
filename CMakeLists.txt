PROJECT(feMorph)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

AUX_SOURCE_DIRECTORY(. SRC_LIST)

# Get Trilinos as one entity
FIND_PACKAGE(Trilinos)

SET( SCH_FLAGS "-Wall -g -O0 -std=c++0x" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SCH_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${SCH_FLAGS}" )

IF(APPLE)
  MESSAGE("APPLE Mode")
  INCLUDE_DIRECTORIES("." "/Users/danschi/Documents/03_software/12_CSparse/include")
  LINK_DIRECTORIES("/Users/danschi/Documents/03_software/12_CSparse/lib")
ELSE()
  MESSAGE("UNIX Mode")
  INCLUDE_DIRECTORIES("." "/home/danschi/Documents/23_Libraries/SuiteSparse/CSparse/Include")
  LINK_DIRECTORIES("/home/danschi/Documents/23_Libraries/SuiteSparse/CSparse/Lib")
ENDIF()

LINK_DIRECTORIES("/opt/local/lib/mpich-gcc46/")

ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_LIST})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} armadillo csparse mpich ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})


