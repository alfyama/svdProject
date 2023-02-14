#!/bin/bash

methods=("dqds" "gr")
types=("float" "double" "long double")

for method in "${methods[@]}"; do
  for type in "${types[@]}"; do
    for (( i=5; i<=100; i+=5 )); do
      filename="test${i}x${i}matrix.csv"
      echo ./svd -method="$method" -type="$type" test/$filename
      ./svd -method="$method" -type="$type" test/$filename
    done
  done
done

for method in "${methods[@]}"; do
  for type in "${types[@]}"; do
    echo ./svd -method="$method" -type="$type" test/test0matrix.csv
    ./svd -method="$method" -type="$type" test/test0matrix.csv
    echo ./svd -method="$method" -type="$type" test/test1matrix.csv
    ./svd -method="$method" -type="$type" test/test1matrix.csv
  done
done
