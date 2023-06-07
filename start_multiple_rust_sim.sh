#!/bin/bash

atoms_numbers_start=400
atoms_numbers_end=402

for var in $(seq $atoms_numbers_start $atoms_numbers_end); do
  echo $var
  mkdir -p ./$var 
  cp ./home_rust_cluster_run.sh ./$var/home_rust_cluster_run.sh
  cd ./$var
  # home_rust_cluster_run.sh $var
  cd ../
done
