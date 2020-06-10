#!/bin/bash

xdir=$(dirname $0)
PSSM_dir=./

if [ $# -le 1 ]; then echo "usage: $0  PSSM_dir *.spd3"; exit 1; fi

PSSM_dir=$1; shift

$xdir/pred_HSE.py $xdir/hsa_full $PSSM_dir NULL $* -hsa
$xdir/pred_HSE.py $xdir/hsb_full $PSSM_dir NULL $* -hsb
