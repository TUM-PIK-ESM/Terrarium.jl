#!/bin/bash
#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --job-name=terrarium-benchmark
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=00:30:00

julia --project=examples --threads=1 test/benchmarks/soil_heat_hydrology_global.jl --device=cpu --samples=10
julia --project=examples --threads=64 test/benchmarks/soil_heat_hydrology_global.jl --device=cpu --samples=10
julia --project=examples --threads=1 test/benchmarks/soil_heat_hydrology_global.jl --device=gpu --samples=10
