# Include Source in Current Folder
AUX_SOURCE_DIRECTORY(. SRC_LIST)

# Use the Boost Libraries
FIND_PACKAGE(Boost 1.40 REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

# Get Trilinos as one entity
IF(useTrilinos)
  FIND_PACKAGE(Trilinos)
  # Include Trilinos Classes
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/trilinos)
ENDIF()

# Add Compilation Flags
SET( SCH_FLAGS "-O3" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SCH_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${SCH_FLAGS}" )

# Add CSparse Include and Lib Folder
IF(useCSparse)
  INCLUDE_DIRECTORIES("." "/home/danschi/Documents/23_Libraries/SuiteSparse/CSparse/Include")
  LINK_DIRECTORIES("/home/danschi/Documents/23_Libraries/SuiteSparse/CSparse/Lib")
ENDIF()

ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_LIST})

IF(useMPI)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} mpich)
ENDIF()

IF(useCSparse)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} csparse)
ENDIF()

IF(useTrilinos)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
ENDIF()

IF(useArmadillo)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} armadillo)
ENDIF()

# Always Include the Boost Library
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES})
