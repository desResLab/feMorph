AUX_SOURCE_DIRECTORY(. SRC_LIST)

# Get Trilinos as one entity
IF(useTrilinos)
  FIND_PACKAGE(Trilinos)
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
