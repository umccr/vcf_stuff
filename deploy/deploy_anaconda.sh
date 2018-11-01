#!/bin/bash
# this script uses the ANACONDA_TOKEN env var. 
# to create a token:
# >>> anaconda login
# >>> anaconda auth -c -n travis --max-age 307584000 --url https://anaconda.org/USERNAME/PACKAGE_NAME --scopes "api:write api:read"
set -e
set -x

PACKAGE_NAME=$1
ANACONDA_TOKEN=$2

if [ -z $PACKAGE_NAME ] ; then
    echo 'PACKAGE_NAME environment variable must be specified' >&2
    exit 1
fi

if [ -z $ANACONDA_TOKEN ] ; then
    echo 'ANACONDA_TOKEN environment variable must be specified' >&2
    exit 1
fi

echo "Converting for macOS..."
conda convert --platform osx-64 $HOME/miniconda/conda-bld/linux-64/${PACKAGE_NAME}-*.tar.bz2 --output-dir $HOME/miniconda/conda-bld/

echo "Deploying to Anaconda.org..."
anaconda -t ${ANACONDA_TOKEN} upload $HOME/miniconda/conda-bld/**/${PACKAGE_NAME}-*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0
