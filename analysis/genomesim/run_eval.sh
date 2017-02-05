#!/bin/sh

if [[ $# != 1 ]]
then
  echo run_eval.sh dir
  exit
fi

DIR=$1

for i in `/bin/ls $DIR/*-out`;
do
  ./eval_err.pl $i > $i.eval
done

grep "#" $DIR/*.eval | sort -nk7
