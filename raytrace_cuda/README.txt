compile raytest.c file int C exec and run

to use cuda file, compile with nvcc + .cu file
(need cuda toolkit installed on your computer to simulate gpu programming)

Notes:
adjust serial to either return lvalues or instantiate vecotrs as lvalue instead of pointers

raytest.c == serial implementation of ray tracking algorithm. Included png of sphere at different completion for correctness

raytrace.cu == cuda implementation of ray tracing algorithm
//  need to include # of blocks on line 195 to test. LEft blank for testing and modifying for trials

10^7 rays && 256x256 grid
2       4       8       16      32      64      128

0.2753  0.2357


serial scaling problem up  on 256x256 sandyb
10^6   2*10^6  4*10^6      8*10^6      10^7

2.385    4.18    8.3376      16.56       20.723

cuda scaling problem up  on 256x256 midway gres=gpu:1
10^6   2*10^6  4*10^6      8*10^6      10^7




(using 16 threads task in sbatch script)
blocksize=1 number of blocks = 1410065408
Cuda,1410065408, 0.221421

blocksize=20 number of blocks = 70503270
Cuda 0.213495

blocksize=40 number of blocks = 35251635
Cuda 0.225487

blocksize=2 number of blocks = 352516352
Cuda 0.249807

blocksize=64 number of blocks = 1101613
Cuda,1410065408, 0.245249

blocksize=256 number of blocks = 110161
Cuda,1410065408, 0.250265

blocksize=2000 number of blocks = 14100
Cuda,1410065408, 0.250653