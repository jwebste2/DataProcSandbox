#!/bin/bash

#
# Set up ROOT paths for Argonne cluster / interactive nodes,
# but only bother doing this if it looks like we are logged into
# an Argonne machine.
#
if [ ! "`echo ${HOSTNAME} | grep anl`" == "" ]; then

    # Required for running ROOT
    echo "BASH : Defining environment variables necessary for ROOT setup"
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    export ALRB_localConfigDir=$HOME/localConfig
    echo "BASH : Setting up cvmfs ROOT"
    source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
    source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh --rootVer=6.04.14-x86_64-slc6-gcc49-opt --skipConfirm

    # Required for processing ProMC files
    echo "BASH : Adding ProMC paths"
    export PROMC=/share/sl6/promc
    export PATH=$PROMC/bin:$PATH
    export PATH=/share/sl6/cmake/bin/:$PATH

fi

#
# Set path to Delphes
#
if [ ! "`echo ${HOSTNAME} | grep anl`" == "" ]; then
    export DELPHESPATH=/users/jwebster/Delphes-3.3.0
else
    export DELPHESPATH=/Users/moon/work/Delphes-3.3.0
fi    
echo "BASH : Setting Delphes path to ${DELPHESPATH}. Please edit the setup script if this is incorrect!"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DELPHESPATH}

#
# define some terminal font styles
# (these will only show up correctly when running interactively, but that's okay)
#
#
echo "BASH : Defining some terminal font styles"
export FONTNONE='\033[00m'
export FONTRED='\033[01;31m'
export FONTGREEN='\033[01;32m'
export FONTYELLOW='\033[01;33m'
export FONTPURPLE='\033[01;35m'
export FONTCYAN='\033[01;36m'
export FONTWHITE='\033[01;37m'
export FONTBOLD='\033[1m'
export FONTUNDERLINE='\033[4m'

echo "BASH : done."
