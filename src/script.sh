#!/bin/bash
make

# create a new log file
date=`date '+%Y_%m_%d_%H_%M_%S'`
file="../timing/test_$date.txt"
echo -e "$file\n"

echo "number_of_light_particles,number_of_medium_particles,number_of_heavy_particles,number_of_particles,number_of_processors,min_time_of_substeps,max_time_of_substeps,avg_time_of_substeps">> $file

# run through timing
subSteps=10
for numParticles in 2000 4000 8000 16000 32000
do
  for p in 32 16 8 4 2
  do
    if (($numParticles % 3))
    then numParticlesLight=$((numParticles/3 + 1))
    else numParticlesLight=$((numParticles/3))
    fi
    if ((($numParticles % 3) == 2))
    then numParticlesHeavy=$((numParticles/3 + 1))
    else numParticlesHeavy=$((numParticles/3))
    fi
    numParticlesMedium=$((numParticles/3))

    echo "mpirun -np $p ./project.x $numParticlesLight $numParticlesMedium $numParticlesHeavy 1 $subSteps 1 1024 1024 \"../images/output\""
    # mpirun -np $p ./project.x $numParticlesLight $numParticlesMedium $numParticlesHeavy 1 $subSteps 1 1024 1024 "../images/output" >> $file
  done
done
