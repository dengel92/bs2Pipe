#!/bin/bash -l

home=$1
models=$2
date=$3

qsub calcmodel.qsub $home $models $date


