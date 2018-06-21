#!/bin/bash

cd ..
DT=`date +%Y%m%d_%H%M%S `
aptly snapshot create t-5_$DT from repo t-5
aptly publish drop xenial
aptly publish snapshot t-5_$DT
rsync -r --progress /home/jh/.aptly/public/ /home/jh/t-5.eu/debian-repo/
