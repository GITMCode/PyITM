#!/bin/bash

PWD=$(pwd)

if [ -e pyitm_GDRIVE ] ; then
  echo "Setting pyitm data path to: $PWD"
  sed "s|PYITMDIR|$PWD|g" pyitm_GDRIVE > pyitm_satmap.csv
  echo "Created pyitm_satmap.csv"
else
  echo "ERROR!! Could not find source satmap file 'pyitm_GDRIVE'"
  echo "Please run this script from PyITM/data/gdrive!"
  exit 1
fi


if [ -e test_GDRIVE ] ; then
  echo "Setting pyitm data path to: $PWD"
  sed "s|PYITMDIR|$PWD|g" test_GDRIVE > test.txt
  echo "Created test.txt"
else
  echo "ERROR!! Could not find source satmap file 'test_GDRIVE'"
  echo "Please run this script from PyITM/data/gdrive!"
  exit 1
fi
    

