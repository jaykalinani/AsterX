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

if [ -z "${SILO_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable SILO_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${SILO_INC_DIRS="${SILO_DIR}/include"}
: ${SILO_LIB_DIRS="${SILO_DIR}/lib"}
: ${SILO_LIBS="siloh5"}

SILO_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${SILO_INC_DIRS})"
SILO_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${SILO_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "SILO_DIR      = ${SILO_DIR}"
echo "SILO_INC_DIRS = ${SILO_INC_DIRS}"
echo "SILO_LIB_DIRS = ${SILO_LIB_DIRS}"
echo "SILO_LIBS     = ${SILO_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(SILO_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(SILO_LIB_DIRS)'
echo 'LIBRARY           $(SILO_LIBS)'
