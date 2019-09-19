set(FIND_QUANTUMTOOLBOX_PATHS
        ../QuantumToolbox)
find_path(QUANTUMTOOLBOX_INCLUDE_DIR Grid.h
        PATHS_SUFFIXES include
        PATHS ${FIND_QUANTUMTOOLBOX_PATHS})
find_library(QUANTUMTOOLBOX_LIBRARY
        NAMES QuantumToolbox
        PATH_SUFFIXES lib
        PATHS ${FIND_QUANTUMTOOLBOX_PATHS})