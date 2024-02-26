#!/usr/bin/env bash

DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUTDIR="./results/${DT}_$1/out"
DIR="./results/${DT}_$1/"

echo "created: $DIR"

mkdir -p $OUTDIR

if [ -d ./out ]; then
  cp -r ./out "$DIR/"
  rm -r ./out/*
fi

if [ -f CRASH ]; then
  mv CRASH "$DIR/out/CRASH"
fi
