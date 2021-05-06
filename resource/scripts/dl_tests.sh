#!/bin/bash

gURL=1gI7_U5DuoNSRKe03PHRrLy5MzLaD7KFU

# match more than 26 word characters  
ggID=$(echo "$gURL" | egrep -o '(\w|-){26,}')

ggURL='https://drive.google.com/uc?export=download'

apt-get install curl
curl -sc /tmp/gcokie "${ggURL}&id=${ggID}" >/dev/null  
getcode="$(awk '/_warning_/ {print $NF}' /tmp/gcokie)"

cmd='curl --insecure -C - -LOJb /tmp/gcokie "${ggURL}&confirm=${getcode}&id=${ggID}"'
echo -e "Downloading from "$gURL"...\n"
cd /MeshTracker/resource
eval $cmd
unzip tests.zip && rm tests.zip
