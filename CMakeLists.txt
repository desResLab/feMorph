# Set Up Project
PROJECT(feMorph)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# ADD OPTIONS
OPTION(useMPI "useMPI" OFF)
OPTION(useCSparse "useCSparse" ON)
OPTION(useTrilinos "useTrilinos" OFF)
OPTION(useArmadillo "useArmadillo" OFF)

# SET MUTUALLY EXCLUSIVE SOLVERS
IF(useTrilinos)
  SET(useArmadillo OFF)
  SET(useCSparse OFF)
ENDIF()
IF(useCSparse)
  SET(useMPI OFF)
  SET(useArmadillo OFF)
  SET(useTrilinos OFF)
ENDIF()
IF(useArmadillo)
  SET(useMPI OFF)
  SET(useCSparse OFF)
  SET(useTrilinos OFF)
ENDIF()

MESSAGE(STATUS "useMPI=${useMPI}")
MESSAGE(STATUS "useCSparse=${useCSparse}")
MESSAGE(STATUS "useTrilinos=${useTrilinos}")
MESSAGE(STATUS "useArmadillo=${useArmadillo}")

# SET DEFINITIONS AND FLAGS
# MPI
if(useMPI)
  ADD_DEFINITIONS("-DUSE_MPI")
endif()
# CSPARSE
if(useCSparse)
  ADD_DEFINITIONS("-DUSE_CSPARSE")
endif()
# TRILINOS
if(useTrilinos)
  ADD_DEFINITIONS("-DUSE_TRILINOS")
endif()
# ARMADILLO
if(useArmadillo)
  ADD_DEFINITIONS("-DUSE_ARMADILLO")
endif()

# ADD SUBDIRS
ADD_SUBDIRECTORY("./src/")
