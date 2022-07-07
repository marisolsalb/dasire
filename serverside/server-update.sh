#!/bin/bash
if ! [ "$HOSTNAME" == "norm.exbio.wzw.tum.de" ]; then exit 1; fi  #This script is designed for norm.

:> /nfs/scratch/amit/DASiRe/serverside-R.tar.gz
pushd /localscratch/marisol/chit/DASiRe/serverside
:> /tmp/dasirebigfiles.txt
find ./R/examples/* -type f >> /tmp/dasirebigfiles.txt
find ./R -name examples -prune -o -name "*.bed.gz" >> /tmp/dasirebigfiles.txt
tar -czvf /nfs/scratch/amit/DASiRe/serverside-R.tar.gz -T /tmp/dasirebigfiles.txt
rsync /nfs/scratch/amit/DASiRe/serverside-R.tar.gz afenn@oscar:~/DASiRe/serverside/serverside-R.tar.gz -av
popd
echo "DASiRe server files are on oscar at ~/DASiRe/serverside/serverside-R.tar.gz"
echo "Go to Oscar and Unzip the file with the commands   ' pushd ~/DASiRe/serverside/ && tar -xvzf serverside-R.tar.gz '"