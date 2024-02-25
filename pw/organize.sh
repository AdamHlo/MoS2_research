#!/usr/bin/env bash

DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUTDIR="./results/$DT/out"
DIR="./results/$DT/"

mkdir -p $OUTDIR

if [ -d ./out ]; then
  cp -r ./out "$DIR/out"
  rm -r ./out/*
fi
