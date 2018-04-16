#!/bin/bash
make

# create a new log file
date=`date '+%Y_%m_%d_%H_%M_%S'`
file="../timing/test_$date.txt"
echo -e "$file\n"

# run through timing
for numParticles in 2000 4000 8000 16000 32000
do
  for p in 32 16 8 4 2 1
  do
    numParticlesLight=$((numParticles/3))
    numParticlesMedium=$((numParticles/3))
    numParticlesHeavy=$((numParticles/3 +1))
    echo "mpirun -np $p ./project.x $numParticlesLight $numParticlesMedium $numParticlesHeavy 1 1 1 1024 1024 \"../images/output\""
    # mpirun -np $p ./project.x $numParticlesLight $numParticlesMedium $numParticlesHeavy 1 1 1 1024 1024 "../images/output" >> $file
  done
done
