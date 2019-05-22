#!/bin/bash

set -e

# trigger build for kkmann/adoptrValidation
body='{
  "request": {
  "branch":"master"
}}'

curl -s -X POST \
   -H "Content-Type: application/json" \
   -H "Accept: application/json" \
   -H "Travis-API-Version: 3" \
   -H "Authorization: token $travis_token" \
   -d "$body" \
   https://api.travis-ci.org/repo/kkmann%2FadoptrValidation/requests

Rscript -e 'pkgdown::deploy_site_github()' &

# Output to the screen every minute to prevent a travis timeout
export PID=$!
while [[ `ps -p $PID | tail -n +2` ]]; do
  echo 'still deploying...'
  sleep 60
done
