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

if [ -z "${YAML_CPP_DIR}" ]; then
    echo "BEGIN ERROR"
    echo "Configuration variable YAML_CPP_DIR is not set"
    echo "END ERROR"
    exit 1
fi

# Set options
: ${YAML_CPP_INC_DIRS="${YAML_CPP_DIR}/include"}
: ${YAML_CPP_LIB_DIRS="${YAML_CPP_DIR}/lib"}
: ${YAML_CPP_LIBS="yaml-cpp"}

YAML_CPP_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${YAML_CPP_INC_DIRS})"
YAML_CPP_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${YAML_CPP_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "YAML_CPP_DIR      = ${YAML_CPP_DIR}"
echo "YAML_CPP_INC_DIRS = ${YAML_CPP_INC_DIRS}"
echo "YAML_CPP_LIB_DIRS = ${YAML_CPP_LIB_DIRS}"
echo "YAML_CPP_LIBS     = ${YAML_CPP_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(YAML_CPP_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(YAML_CPP_LIB_DIRS)'
echo 'LIBRARY           $(YAML_CPP_LIBS)'
