#!/bin/bash
################################################################################
function script () {
cat <<EOF
#!/bin/bash
#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM -L vnode=1
#PJM -L vnode-core=1
#PJM -L elapse=1h

python ../EMh.py --cmf=$cmf --param="$param" --a="$a" --b="$b" --method=$method --Mmax=$Mmax --hmin=$hmin --hmax=$hmax --ymin=$ymin --prec=$prec
EOF
}
################################################################################
function mkdir_run () {
dir=$cmf$param\_a$a\_b$b\_EMh_$method\_prec$prec
if [ -d $dir ]; then
  rm -rf $dir
fi

mkdir $dir
cd $dir
echo $dir
script >> run.sh
pjsub run.sh
cd ..
}
################################################################################
cmf=inverse_power
param_list=("0.5" "1" "2")
b=1
method="qf"
Mmax=17

a="0.5"; hmin="0.05"; hmax="20"; ymin="1e-50"; prec_list=(248 500)
for param in "${param_list[@]}"; do
  for prec in "${prec_list[@]}"; do
    mkdir_run
  done
done

a="0.0009765625"; hmin="5e-2"; hmax="2e4"; ymin="1e-30"; prec_list=(184 500)
for param in "${param_list[@]}"; do
  for prec in "${prec_list[@]}"; do
    mkdir_run
  done
done

