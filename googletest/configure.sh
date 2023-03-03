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

# Set options
GOOGLETEST_DIR=$(dirname "$0")
: ${GOOGLETEST_INC_DIRS="$GOOGLETEST_DIR/include"}
: ${GOOGLETEST_LIB_DIRS=""}
: ${GOOGLETEST_LIBS=''}

GOOGLETEST_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${GOOGLETEST_INC_DIRS})"
GOOGLETEST_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${GOOGLETEST_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "GOOGLETEST_DIR      = ${GOOGLETEST_DIR}"
echo "GOOGLETEST_INC_DIRS = ${GOOGLETEST_INC_DIRS}"
echo "GOOGLETEST_LIB_DIRS = ${GOOGLETEST_LIB_DIRS}"
echo "GOOGLETEST_LIBS     = ${GOOGLETEST_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(GOOGLETEST_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(GOOGLETEST_LIB_DIRS)'
echo 'LIBRARY           $(GOOGLETEST_LIBS)'
