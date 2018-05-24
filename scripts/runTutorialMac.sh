#!/bin/bash

# XQuartz preferences -> Security tab -> Turn on 'Allow connection from network clients'
# Restart XQuartz

docker pull trausch/variant-calling
IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
docker run -d -it -e DISPLAY=${IP}:0 -u ubuntu -v /tmp/.X11-unix/:/tmp/.X11-unix:ro trausch/variant-calling /bin/bash
export containerId=`docker ps -l -q`
xhost + ${IP}
echo "Press Enter to start container..."
docker start -i $containerId
