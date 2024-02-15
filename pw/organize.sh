#!/usr/bin/env bash

DT=`date +"%d.%b.%Y_%T"`
OUTDIR="./results/$DT/out"
DIR="./results/$DT/"

mkdir -p $OUTDIR

if [ -f ./electron_minimization.out ] && [ -f ./relaxation.out ]; then
  echo "Both electron_minimization.out and relaxation.out found, please check if you have proper setup."
  exit 1
fi

if [ -f ./electron_minimization.out ]; then
   mv ./electron_minimization.out "$DIR/electron_minimization.out"
fi

if [ -f ./relaxation.out ]; then
   mv ./relaxation.out "$DIR/relaxation.out"
fi

if [ -f ./out/MoS2.xml ]; then
   mv ./out/MoS2.xml "$OUTDIR/MoS2.xml"
fi

if [ -d ./out/MoS2.save ]; then
   mv ./out/MoS2.save "$OUTDIR/MoS2.save"
fi
