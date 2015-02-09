#!/bin/sh
#PBS -N Feynman_Simulator.MC
#PBS -o /home/chen-huang-zhang/Feynman_Simulator/Output
#PBS -e /home/chen-huang-zhang/Feynman_Simulator/Error
cd /home/chen-huang-zhang/Feynman_Simulator
echo $PBS_JOBID >> id_job.log
/home/chen-huang-zhang/Feynman_Simulator/simulator.exe -f /home/chen-huang-zhang/Feynman_Simulator/infile/_in_MC_0
