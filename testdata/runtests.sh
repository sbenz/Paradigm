#!/bin/bash

echo Testing node splitting, should take less than a minute
./pathwaytab2daifg needs_split.pathway.tab needs_split.cfg  \
    | diff needs_split.out - || exit 1

echo Testing inference, should take less than a minute
../hgFactorGraph -c noem.cfg -p small_pid_66_pathway.tab -b small_pid_66 \
    | python ../helperScripts/diffSwarmFiles.py noem.cfg.out -\
    | diff - /dev/null \
    || exit 1

echo Testing EM, should take approximately five minutes
../hgFactorGraph -c em_simple.cfg -p small_pid_66_pathway.tab -b small_pid_66 \
    python ../helperScripts/diffSwarmFiles.py em_simple.cfg.out -\
    diff - /dev/null || exit 1

