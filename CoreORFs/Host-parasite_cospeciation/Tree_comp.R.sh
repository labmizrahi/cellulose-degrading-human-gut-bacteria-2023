#!/bin/bash
#$ -S /bin/bash


Rscript /cellulose-degrading-human-gut-bacteria-2023/CoreORFs/Host-parasite_cospeciation/Tree_comp.R \
    --host_tree '/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/Host-parasite_cospeciation/Host.nwk' \
    --input_dir "/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/197_Trees_Rooted_with_OutGroup/" \
    --bac2host_file "/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/Host-parasite_cospeciation/Bac2Host" \
    --tree_file_pattern ".treefile" \
    --lable2drop 'Clostridium_thermobutyricum' \
    --tree_outgroup 'Clostridium_thermobutyricum' \
    --nboot 1000 \
    --method "mantel" \
    --cpu 1 \
    --output_dir '/' 
