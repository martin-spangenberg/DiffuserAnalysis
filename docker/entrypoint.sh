#!/bin/bash

source /root/bin/thisroot.sh

export PYTHONPATH=/:$PYTHONPATH
export PYTHONPATH=/plots:$PYTHONPATH

python3 -m "scripts.${1}" "${@:2}"

#cp *.png /plots/

#python3 "$@"

#echo "$@"
#echo "${@:2}"
