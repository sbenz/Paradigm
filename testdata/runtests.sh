#!/bin/bash

set -o pipefail

echo Testing node splitting [1/2], should take less than a minute
../pathwaytab2daifg needs_split_1.pathway.tab needs_split_1.cfg  \
    | diff needs_split_1.out - || exit 1

echo Testing node splitting [2/2], should take less than a minute
../pathwaytab2daifg needs_split_2.pathway.tab needs_split_2.cfg  \
    | diff needs_split_2.out - || exit 1

echo Testing inference, should take less than a minute
/usr/bin/time ../hgFactorGraph -c noem.cfg -p small_pid_66_pathway.tab -b small_pid_66 \
    | python ../helperScripts/diffSwarmFiles.py noem.cfg.out -\
    | diff - /dev/null \
    || exit 1

echo Testing EM, should take approximately five minutes
/usr/bin/time ../hgFactorGraph -c em_simple.cfg -p small_pid_66_pathway.tab -b small_pid_66 \
    | python ../helperScripts/diffSwarmFiles.py em_simple.cfg.out -\
    | diff - /dev/null \
    || exit 1
