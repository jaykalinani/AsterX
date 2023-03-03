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

if [ -z "${OPENPMD_API_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable OPENPMD_API_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${OPENPMD_API_INC_DIRS="${OPENPMD_API_DIR}/include"}
: ${OPENPMD_API_LIB_DIRS="${OPENPMD_API_DIR}/lib"}
: ${OPENPMD_API_LIBS="openPMD"}

OPENPMD_API_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${OPENPMD_API_INC_DIRS})"
OPENPMD_API_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${OPENPMD_API_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "OPENPMD_API_DIR      = ${OPENPMD_API_DIR}"
echo "OPENPMD_API_INC_DIRS = ${OPENPMD_API_INC_DIRS}"
echo "OPENPMD_API_LIB_DIRS = ${OPENPMD_API_LIB_DIRS}"
echo "OPENPMD_API_LIBS     = ${OPENPMD_API_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(OPENPMD_API_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(OPENPMD_API_LIB_DIRS)'
echo 'LIBRARY           $(OPENPMD_API_LIBS)'
