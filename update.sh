#!/bin/bash

git add $1
git commit -m \"$2\"
git push

rsync -Pavz $1 sourav.c@www.tifr.res.in:/home/depts/sourav.c/public_html/
