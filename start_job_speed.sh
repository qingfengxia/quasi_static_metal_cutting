ORIGIN="parameter_template.py"
SCRIPT="parameter.py"
VAR="cutting_speed"
#first of all, comment out the var assign if after anchor line (##ZZ##) in parameter_template
#also check result file name

#double quote to expand variable into string
for VALUE in 0.1 0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0; do
#for VALUE in 2.5 3.0 3.5 4.0; do
    newline="$VAR=$VALUE"
    echo "processing parameter:$newline"
    sed "s/##ZZ##/$newline/" $ORIGIN > $SCRIPT
    #is_batch_mode = False
    python3 metal_cut_ht_6.py
done



