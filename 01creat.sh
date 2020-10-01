#!/bin/bash

dirname=1000_frame
j=0
for i in {1..800}
do
  j=$((j+1))
  echo $i $j
  ./a.out  $dirname/water_$i.xyz 
  mv $dirname/water_$i.xyz.ft Train_data/$j.ft
  mv $dirname/water_$i.xyz.wc Train_data/$j.wc
  rm $dirname/water_$i.xyz_labframe.xyz $dirname/water_$i.xyz_molframe.xyz $dirname/water_$i.xyz_nbratoms.xyz
done 

j=0
for i in {801..999}
do
  j=$((j+1))
  echo $i $j
  ./a.out  $dirname/water_$i.xyz
  mv $dirname/water_$i.xyz.ft Test_data/$j.ft
  mv $dirname/water_$i.xyz.wc Test_data/$j.wc
  rm $dirname/water_$i.xyz_labframe.xyz $dirname/water_$i.xyz_molframe.xyz $dirname/water_$i.xyz_nbratoms.xyz
done 
