#!/bin/bash

 mpirun -np 1 --hostfile hostfile.txt main 50000 20 >> test50kX20np1.txt
echo "Test 1 Terminato"
 mpirun -np 2 --hostfile hostfile.txt main 50000 20 >> test50kX20np2.txt
echo "Test 2 terminato"
 mpirun -np 4 --hostfile hostfile.txt main 50000 20 >> test50kX20np4.txt
echo "Test 4 terminato"
 mpirun -np 6 --hostfile hostfile.txt main 50000 20 >> test50kX20np6.txt
echo "Test 6 terminato"
 mpirun -np 8 --hostfile hostfile.txt main 50000 20 >> test50kX20np8.txt
echo "Test 8 terminato"
mpirun -np 10 --hostfile hostfile.txt main 50000 20 >> test50kX20np10.txt
echo "Test 10 terminato"
 mpirun -np 12 --hostfile hostfile.txt main 50000 20 >> test50kX20np12.txt
echo "Test 12 terminato"
 mpirun -np 14 --hostfile hostfile.txt main 50000 20 >> test50kX20np14.txt
echo "Test 14 terminato"
 mpirun -np 16 --hostfile hostfile.txt main 50000 20 >> test50kX20np16.txt
echo "Test 16 terminato"