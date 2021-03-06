#!/bin/bash
## build_containers.sh 
## Singularity images of HPV pipeline
##
## Copyright (c) 2019-2020 Institut Curie
## Author(s): Nicolas Servant, Philippe La Rosa
## Contact: nicolas.servant@curie.fr, philippe.larosa@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details
##
##
## Run by root user or user with sudo
##
## usage : bash build_containers.sh 2>&1 | tee -a build_containers.log 
##

mkdir -p images
for RECI in $(ls recipes/*.def)
 do 
   IMGNAME=$(basename ${RECI} .def)
   echo "## build image ${IMGNAME}" 
   sudo singularity build images/${IMGNAME}.simg ${RECI}
 done 

