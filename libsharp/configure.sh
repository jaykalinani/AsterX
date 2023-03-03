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

if [ -z "${LIBSHARP_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable LIBSHARP_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${LIBSHARP_INC_DIRS="${LIBSHARP_DIR}/include"}
: ${LIBSHARP_LIB_DIRS="${LIBSHARP_DIR}/lib"}
: ${LIBSHARP_LIBS="sharp fftpack c_utils"}

LIBSHARP_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${LIBSHARP_INC_DIRS})"
LIBSHARP_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${LIBSHARP_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "LIBSHARP_DIR      = ${LIBSHARP_DIR}"
echo "LIBSHARP_INC_DIRS = ${LIBSHARP_INC_DIRS}"
echo "LIBSHARP_LIB_DIRS = ${LIBSHARP_LIB_DIRS}"
echo "LIBSHARP_LIBS     = ${LIBSHARP_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(LIBSHARP_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(LIBSHARP_LIB_DIRS)'
echo 'LIBRARY           $(LIBSHARP_LIBS)'
