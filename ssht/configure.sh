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

if [ -z "${SSHT_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable SSHT_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${SSHT_INC_DIRS="${SSHT_DIR}/include"}
: ${SSHT_LIB_DIRS="${SSHT_DIR}/lib"}
: ${SSHT_LIBS="ssht"}

SSHT_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${SSHT_INC_DIRS})"
SSHT_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${SSHT_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "SSHT_DIR      = ${SSHT_DIR}"
echo "SSHT_INC_DIRS = ${SSHT_INC_DIRS}"
echo "SSHT_LIB_DIRS = ${SSHT_LIB_DIRS}"
echo "SSHT_LIBS     = ${SSHT_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(SSHT_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(SSHT_LIB_DIRS)'
echo 'LIBRARY           $(SSHT_LIBS)'
