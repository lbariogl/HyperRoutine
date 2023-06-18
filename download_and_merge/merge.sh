#!/bin/bash

INPUT=$1
OUTPUT=$2
ISMC=${3:-false}


root -l -b << EOF
.L MergeTrees.cc++
MergeTrees("${INPUT}", "${OUTPUT}", ${ISMC})
EOF
