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

if [ -z "${NSIMD_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable NSIMD_DIR is not set"
    echo "END ERROR"
    exit 1
fi

if [ -z "${NSIMD_ARCH}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable NSIMD_ARCH is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${NSIMD_INC_DIRS="${NSIMD_DIR}/include"}
: ${NSIMD_LIB_DIRS="${NSIMD_DIR}/lib"}
: ${NSIMD_LIBS="nsimd_${NSIMD_ARCH}"}

NSIMD_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${NSIMD_INC_DIRS})"
NSIMD_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${NSIMD_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN DEFINE"
for opt in ${NSIMD_OPTIONS}; do
    echo "NSIMD_${opt}"
done
echo "END DEFINE"

echo "BEGIN MAKE_DEFINITION"
echo "NSIMD_DIR      = ${NSIMD_DIR}"
echo "NSIMD_INC_DIRS = ${NSIMD_INC_DIRS}"
echo "NSIMD_LIB_DIRS = ${NSIMD_LIB_DIRS}"
echo "NSIMD_LIBS     = ${NSIMD_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(NSIMD_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(NSIMD_LIB_DIRS)'
echo 'LIBRARY           $(NSIMD_LIBS)'
