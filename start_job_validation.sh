#!/bin/bash
ORIGIN="parameter_template.py"
SCRIPT="parameter.py"
VAR="case_id"
#first of all, comment out the case_id in parameter_template

#double quote to expand variable into string
for VALUE in 1 2 3 4 5 6 7 8; do
    newline="$VAR=$VALUE"
    echo "processing parameter:$newline"
    sed "s/##ZZ##/$newline/" $ORIGIN > $SCRIPT
    #is_batch_mode = False
    python3 metal_cut_ht_6.py
done

