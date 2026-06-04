#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build
cmake --build build -j8
mkdir -p out
./build/ogden_3d_uniaxial > out/uniaxial.csv
./build/ogden_surface_triangle > out/surface_triangle.txt

echo "Wrote out/uniaxial.csv"
echo "Wrote out/surface_triangle.txt"
