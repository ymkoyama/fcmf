#!/bin/bash
################################################################################
function print_error () {
prefix=$cmf$param\_a$a\_b$b\_EMh_$method
dir0=$prefix\_prec$prec_ref
dir1=$prefix\_prec$prec
echo "target    $dir1"
echo "reference $dir0"

for M in `seq $Mmax`; do
  name=`printf "EMh_f0_M%04d.dat" $M`
  file0=$dir0/$name
  file1=$dir1/$name
  python print_error.py --file0=$file0 --file1=$file1 --prec=$prec_ref | awk -v M=$M '$1 == "max_rel_error"{printf("M %4d MRE(h) %s MRE(EMh/f0) %s\n", M, $2, $3)}'
done
}
################################################################################
cmf=inverse_power
param_list=("0.5" "1" "2")
b=1
method="qf"
Mmax=17
prec_ref=500

a="0.5"; prec=248
for param in "${param_list[@]}"; do
  print_error
done

a="0.0009765625"; prec=184
for param in "${param_list[@]}"; do
  print_error
done

