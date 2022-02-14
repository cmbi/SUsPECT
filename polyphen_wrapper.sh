#!/bin/bash
M=15   # number of program instances to run
for (( N=1; N<=$M; N++ )); do
  /var/lib/polyphen-2.2.2/bin/run_pph.pl -r $N/$M "$@" 1>pph$N.features 2>pph$N.log &
done
wait
rm -f pph.features pph.log
for (( N=1; N<=$M; N++ )); do
  cat pph$N.features >>pph.features
  cat pph$N.log >>pph.log
  rm -f pph$N.features pph$N.log
done