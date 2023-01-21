#!/bin/bash
################################################################################
function mkdir_run () {
dir=$name
if [ -d $dir ]; then
  rm -rf $dir
fi

mkdir $dir
cd $dir
echo $dir

for k in "${!MDS_list[@]}"; do
  method=${method_list[k]}
  MDS=${MDS_list[k]}
  phi=${phi_list[k]}
  prefix=../$name\_$method\_MDS$MDS\_prec$prec/$method\_$phi
  echo $prefix >> prefix.dat
done

for M in `seq 1 $Mmax`; do
  python -u ../fig_EMx.py --prefix=prefix.dat --M=$M | tee -a fig_EMx.log
done

out=$dir\_EMx.pdf
pdfunite E_M*.pdf $out
echo $dir/$out

cd ..
}
################################################################################
method_list=("best" "gauss" "gauss" "gauss" "gauss" "gauss")
phi_list=("P2" "Phi" "exp" "P2" "P1" "R0_1")

prec=248; MDS_list=(96 96 96 96 96 96); Mmax=17
name=inverse_power0.5_a0.5_b1; mkdir_run
name=inverse_power1_a0.5_b1; mkdir_run
name=inverse_power2_a0.5_b1; mkdir_run

prec=184; MDS_list=(192 1536 1536 1536 1536 1536); Mmax=17
name=inverse_power0.5_a0.0009765625_b1; mkdir_run
name=inverse_power1_a0.0009765625_b1; mkdir_run
name=inverse_power2_a0.0009765625_b1; mkdir_run

