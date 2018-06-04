#!/bin/bash

 mpirun -np 1 --hostfile hostfile.txt main 10000 20 >> weakT1_10kX20np1.txt
echo "Weak Test 1: -np 1 Terminato"
 mpirun -np 2 --hostfile hostfile.txt main 20000 20 >> weakT1_10kX20np2.txt
echo "Weak Test 1: -np 2 terminato"
 mpirun -np 4 --hostfile hostfile.txt main 40000 20 >> weakT1_10kX20np4.txt
echo "Weak Test 1: -np 4 terminato"
 mpirun -np 6 --hostfile hostfile.txt main 60000 20 >> weakT1_10kX20np6.txt
echo "Weak Test 1: -np 6 terminato"
 mpirun -np 8 --hostfile hostfile.txt main 80000 20 >> weakT1_10kX20np8.txt
echo "Weak Test 1: -np 8 terminato"
mpirun -np 10 --hostfile hostfile.txt main 100000 20 >> weakT1_10kX20np10.txt
echo "Weak Test 1: -np 10 terminato"
 mpirun -np 12 --hostfile hostfile.txt main 120000 20 >> weakT1_10kX20np12.txt
echo "Weak Test 1: -np 12 terminato"
 mpirun -np 14 --hostfile hostfile.txt main 140000 20 >> weakT1_10kX20np14.txt
echo "Weak Test 1: -np 14 terminato"
 mpirun -np 16 --hostfile hostfile.txt main 160000 20 >> weakT1_10kX20np16.txt
echo "Weak Test 1: -np 16 terminato"
 mpirun -np 1 --hostfile hostfile.txt main 3000 20 >> weakT2_3kX20np1.txt
echo "Weak Test 2: -np 1 Terminato"
 mpirun -np 2 --hostfile hostfile.txt main 6000 20 >> weakT2_3kX20np2.txt
echo "Weak Test 2: -np 2 terminato"
 mpirun -np 4 --hostfile hostfile.txt main 12000 20 >> weakT2_3kX20np4.txt
echo "Weak Test 2: -np 4 terminato"
 mpirun -np 6 --hostfile hostfile.txt main 18000 20 >> weakT2_3kX20np6.txt
echo "Weak Test 2: -np 6 terminato"
 mpirun -np 8 --hostfile hostfile.txt main 24000 20 >> weakT2_3kX20np8.txt
echo "Weak Test 2: -np 8 terminato"
mpirun -np 10 --hostfile hostfile.txt main 30000 20 >> weakT2_3kX20np10.txt
echo "Weak Test 2: -np 10 terminato"
 mpirun -np 12 --hostfile hostfile.txt main 36000 20 >> weakT2_3kX20np12.txt
echo "Weak Test 2: -np 12 terminato"
 mpirun -np 14 --hostfile hostfile.txt main 42000 20 >> weakT2_3kX20np14.txt
echo "Weak Test 2: -np 14 terminato"
 mpirun -np 16 --hostfile hostfile.txt main 48000 20 >> weakT2_3kX20np16.txt
echo "Weak Test 2: -np 16 terminato"