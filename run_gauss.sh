#!/bin/bash
################################################################################
function run () {
cat <<EOF
#!/bin/bash
#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM -L vnode=1
#PJM -L vnode-core=1
#PJM -L elapse=1h

python ../expsum.py --cmf=$cmf --param=$param --a=$a --b=$b --method=$method --Mmax=$Mmax --phi=Phi,exp,P2,P1,R0_1 --MDS=$MDS --prec=$prec
EOF
}
################################################################################
function mkdir_run () {
dir=$cmf$param\_a$a\_b$b\_$method\_MDS$MDS\_prec$prec
if [ -d $dir ]; then
  rm -rf $dir
fi

mkdir $dir
cd $dir
echo $dir
run >> run.sh
pjsub run.sh
cd ..
}
################################################################################
cmf=inverse_power
param_list=("0.5" "1" "2")
b=1
method="gauss"
Mmax=17

a="0.5"; prec=248; MDS_list=(24 48 96 192)
for param in "${param_list[@]}"; do
  for MDS in "${MDS_list[@]}"; do
    mkdir_run
  done
done

a="0.0009765625"; prec=184; MDS_list=(24 48 96 192 384 768 1536 3072)
for param in "${param_list[@]}"; do
  for MDS in "${MDS_list[@]}"; do
    mkdir_run
  done
done

