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

if [ -z "${KOKKOS_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable KOKKOS_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${KOKKOS_INC_DIRS="${KOKKOS_DIR}/include"}
: ${KOKKOS_LIB_DIRS="${KOKKOS_DIR}/lib"}
: ${KOKKOS_LIBS="kokkos"}

KOKKOS_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${KOKKOS_INC_DIRS})"
KOKKOS_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${KOKKOS_LIB_DIRS})"

echo "BEGIN MAKE_DEFINITION"
echo "KOKKOS_DIR      = ${KOKKOS_DIR}"
echo "KOKKOS_INC_DIRS = ${KOKKOS_INC_DIRS}"
echo "KOKKOS_LIB_DIRS = ${KOKKOS_LIB_DIRS}"
echo "KOKKOS_LIBS     = ${KOKKOS_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(KOKKOS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(KOKKOS_LIB_DIRS)'
echo 'LIBRARY           $(KOKKOS_LIBS)'
