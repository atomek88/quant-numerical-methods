HPC MPCS 51087 - Tomasz Michalik
Prob set 4 - Poisson Parallel Solver

# Serial implementation ****(used most of what provided code showed)****
Running on 75x75 physical problem size:
Dense method memory usage; 241.61MB (says for program), yet calculation shows 0.064 MB
Dense method runtime: 50.05 sec
For 10,000 x 10,000 problem, memory =  17592186029767.20 MB

Sparse Memory usage = 0.21MB (says for program), yet calculation shows 0.003
Sparse runtime = 0.07 seconds
For 10,000 x 10,000 problem, memory =  3814.70 MB
serial

512 x 512 matrix

parallel
1:
2 : 9.3sec
4: 4.77sec
8: 2.48 sec
16: 1.37sec  1561 iterations
32: 0.95sec
64: 1.1sec

weak scaling

n=512  2
9.3sec
n=720  4
13.000sec
n=1024 8
19.7sec
n=1440 16
30.4sec
n=2048 32
42.92sec
n=2880 64
67.4sec