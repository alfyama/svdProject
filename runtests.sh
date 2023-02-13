#!/bin/bash

methods=("gr" "dqds")
types=("float" "double" "long double")

for method in "${methods[@]}"; do
  for type in "${types[@]}"; do
    for (( i=5; i<=50; i+=5 )); do
      filename="testmatrix_${i}x${i}.csv"
      ./svd -method="$method" -type="$type" tests/$filename
      echo $filename
    done
  done
done
