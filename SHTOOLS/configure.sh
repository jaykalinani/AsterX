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

if [ -z "${SHTOOLS_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable SHTOOLS_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${SHTOOLS_INC_DIRS="${SHTOOLS_DIR}/include"}
: ${SHTOOLS_LIB_DIRS="${SHTOOLS_DIR}/lib"}
: ${SHTOOLS_LIBS="SHTOOLS-mp"}

SHTOOLS_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${SHTOOLS_INC_DIRS})"
SHTOOLS_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${SHTOOLS_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "SHTOOLS_DIR      = ${SHTOOLS_DIR}"
echo "SHTOOLS_INC_DIRS = ${SHTOOLS_INC_DIRS}"
echo "SHTOOLS_LIB_DIRS = ${SHTOOLS_LIB_DIRS}"
echo "SHTOOLS_LIBS     = ${SHTOOLS_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(SHTOOLS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(SHTOOLS_LIB_DIRS)'
echo 'LIBRARY           $(SHTOOLS_LIBS)'
