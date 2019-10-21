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

if [ -z "${AMREX_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable AMREX_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${AMREX_INC_DIRS="${AMREX_DIR}/include"}
: ${AMREX_LIB_DIRS="${AMREX_DIR}/lib"}
: ${AMREX_LIBS='amrex'}

AMREX_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${AMREX_INC_DIRS})"
AMREX_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${AMREX_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "AMREX_DIR      = ${AMREX_DIR}"
echo "AMREX_INC_DIRS = ${AMREX_INC_DIRS}"
echo "AMREX_LIB_DIRS = ${AMREX_LIB_DIRS}"
echo "AMREX_LIBS     = ${AMREX_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(AMREX_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(AMREX_LIB_DIRS)'
echo 'LIBRARY           $(AMREX_LIBS)'
