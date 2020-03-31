#!/bin/bash

function run_scan() {
  # bash -x scripts/indirect_constraint_contour.sh $1 $2 26
  bash -x scripts/all_constraint_contour.sh $1 $2 26
}


for lam2 in "-1.66 1.16" "-1.16 -0.66" "-0.66 -0.16" "-0.16 0.34" "0.34 0.84" "0.84 1.34"; do
 low=$(echo $lam2 | awk '{print $1}')
 high=$(echo $lam2 | awk '{print $2}')

 run_scan $low $high > log_run_scan_all_${low}_${high}.out 2>&1 & 
 disown
done

