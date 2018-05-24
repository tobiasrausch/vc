#!/bin/bash

docker pull trausch/variant-calling
docker run -d -it -e DISPLAY=${DISPLAY} -u ubuntu -v /tmp/.X11-unix/:/tmp/.X11-unix:ro trausch/variant-calling /bin/bash
export containerId=`docker ps -l -q`
xhost +local:`docker inspect --format='{{ .Config.Hostname }}' $containerId`
echo "Press Enter to start container..."
docker start -i $containerId
