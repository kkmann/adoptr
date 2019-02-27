#!/bin/bash

Rscript -e 'pkgdown::deploy_site_github()' &

# Output to the screen every minute to prevent a travis timeout
export PID=$!
while [[ `ps -p $PID | tail -n +2` ]]; do
  echo 'still deploying...'
  sleep 60
done
