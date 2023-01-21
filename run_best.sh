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

python ../expsum.py --cmf=$cmf --param=$param --a=$a --b=$b --method=$method --Mmin=$M --Mmax=$M --phi=$phi --MDS=$MDS --prec=$prec
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

for M in `seq $Mmax`; do
  runfile=`printf "best%04d.sh" $M`
  echo $runfile
  run >> $runfile
  pjsub $runfile
done

cd ..
}
################################################################################
cmf=inverse_power
param_list=("0.5" "1" "2")
b=1
method=best
phi=P2
Mmax=17

a="0.5"; MDS=96; prec=248
for param in "${param_list[@]}"; do
  mkdir_run
done

a="0.0009765625"; MDS=192; prec=184
for param in "${param_list[@]}"; do
  mkdir_run
done

