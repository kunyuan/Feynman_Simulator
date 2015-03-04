#!/bin/bash
file=status.ipynb
if [ -n "$1" ]
  then
  if [ $1 = "-l" ] || [ $1 = "--local" ]; then
      ipython notebook status.ipynb
  fi
  else
      ipython notebook --no-browser --ip=* --port=8888
      #http://ip-address:8888/notebooks/status.ipynb to access the file
fi
