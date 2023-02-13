#!/bin/bash

methods=("dqds" "gr")
types=("float" "double" "long double")

for method in "${methods[@]}"; do
  for type in "${types[@]}"; do
    for (( i=5; i<=50; i+=5 )); do
      filename="testmatrix_${i}x${i}.csv"
      echo ./svd -method="$method" -type="$type" test/$filename
      ./svd -method="$method" -type="$type" test/$filename
    done
  done
done