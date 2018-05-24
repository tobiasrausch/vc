#!/bin/bash

# Install Docker Toolbox
# Install XMing
# Use XLaunch -> Multiple Windows -> Start no client -> Check "No Access Control"
# Launch Docker Quickstart Terminal

docker pull trausch/variant-calling
IP=$(ipconfig | grep "IPv4 Address" | head -n 1 | sed 's/^.* //')
docker run -d -it -e DISPLAY=${IP}:0 -u ubuntu -v /tmp/.X11-unix/:/tmp/.X11-unix:ro trausch/variant-calling /bin/bash
export containerId=`docker ps -l -q`
echo "Press Enter to start container..."
docker start -i $containerId
