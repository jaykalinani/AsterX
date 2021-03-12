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

if [ -z "${CUDA_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable CUDA_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${CUDA_INC_DIRS="${CUDA_DIR}/include"}
: ${CUDA_LIB_DIRS="${CUDA_DIR}/lib"}
: ${CUDA_LIBS="CUDA"}

CUDA_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${CUDA_INC_DIRS})"
CUDA_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${CUDA_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "CUDA_DIR      = ${CUDA_DIR}"
echo "CUDA_INC_DIRS = ${CUDA_INC_DIRS}"
echo "CUDA_LIB_DIRS = ${CUDA_LIB_DIRS}"
echo "CUDA_LIBS     = ${CUDA_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(CUDA_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(CUDA_LIB_DIRS)'
echo 'LIBRARY           $(CUDA_LIBS)'
