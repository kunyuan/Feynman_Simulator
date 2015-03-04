#!/bin/bash
file=status.ipynb
if [ -n "$1" ]
  then
  if [ $1 = "-l" ] || [ $1 = "--local" ]; then
      ipython notebook status.ipynb
  fi
  else
      ipython notebook --no-browser --ip=* --port=8888 status.ipynb
fi
