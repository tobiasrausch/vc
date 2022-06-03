#!/bin/bash

if [ -d data/sv/ ]
then
    cat README.md | awk '/^```bash$/,/^```$/  {print} {next}' | grep -v '^`' | grep -v "^igv" > run.sh
    chmod a+x run.sh
    ./run.sh
fi
