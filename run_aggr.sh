#!/bin/bash -l

date=$1
mainoutputdir=$2

qsub centaggr.qsub $date $mainoutputdir


