#!/bin/bash

INPUT_DIR=$1
FILE_TYPE=$2

COUNTER=0

for file in $(alien.py find ${INPUT_DIR} ${FILE_TYPE}); do

  if [ file != ${INPUT_DIR}/${FILE_TYPE} ]; then
    alien.py cp alien://${file} file://${FILE_TYPE::${#FILE_TYPE}-5}_${COUNTER}.root
    COUNTER=$((COUNTER++))
  fi
done

hadd ${FILE_TYPE} ${FILE_TYPE::${#FILE_TYPE}-5}_*
rm ${FILE_TYPE::${#FILE_TYPE}-5}_*
