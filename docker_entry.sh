#!/bin/bash
set -e

if [ $# -eq 0 ]
  then
    jupyter lab --ip=0.0.0.0 --NotebookApp.token='local-development' --allow-root --no-browser &> /dev/null &
    code-server --allow-http --no-auth --data-dir /src
  else
    exec "cli/${@}"
fi
