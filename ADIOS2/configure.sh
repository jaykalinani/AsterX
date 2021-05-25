#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Configure Cactus
################################################################################

if [ -z "${ADIOS2_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable ADIOS2_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${ADIOS2_INC_DIRS="${ADIOS2_DIR}/include"}
: ${ADIOS2_LIB_DIRS="${ADIOS2_DIR}/lib"}
: ${ADIOS2_LIBS="adios2_cxx11_mpi adios2_cxx11"}

ADIOS2_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${ADIOS2_INC_DIRS})"
ADIOS2_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${ADIOS2_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN DEFINE"
echo "ADIOS2_USE_MPI 1"
echo "END DEFINE"

echo "BEGIN MAKE_DEFINITION"
echo "ADIOS2_DIR      = ${ADIOS2_DIR}"
echo "ADIOS2_INC_DIRS = ${ADIOS2_INC_DIRS}"
echo "ADIOS2_LIB_DIRS = ${ADIOS2_LIB_DIRS}"
echo "ADIOS2_LIBS     = ${ADIOS2_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(ADIOS2_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(ADIOS2_LIB_DIRS)'
echo 'LIBRARY           $(ADIOS2_LIBS)'
