Tomasz Michalik Problem set 3 - Julia set and latency tests
Notes:
**had issues with ostream compiling on midway. For dynamic arg #2 = 1, For static arg #2 = 0; chunk arg #3 optional for dynamic

** final Julia.cpp final had issues compiling on my computer due to linking error so did not include exe but compiled fine on midway

** included couple sbatch scripts, modified #ranks and static vs dynamic  for scaling tests

** static vs dynamic uses different print routines due to variability in dimensions for each rank of dyna/static

pingtest = used for latency
pingerB = used for bandwidth

avg time on midway bandwidth (intra):
1KB  :  277768 kb/sec   277 mb/sec      0.27 gb/sec
1MB  :  3381276 kb/sec  3381 mb/sec     3.38 gb/sec
1GB  :  3883915 kb/sec  3883 mb/sec     3.88 gb/sec

internode bandwidth
1KB  :  120353 kb/sec   120 mb/sec      0.12 gb/sec
1MB  :  3429240 kb/sec  3429 mb/sec     3.43 gb/sec
1GB  :  3587871 kb/sec  3587 mb/sec     3.58 gb/sec


time on midway latency:
Inter node: avg = 3.0347 avg microSec for 1 int (4 bytes)
Intra node: avg = 1.03 microseconds

possible print routine
// if rank == i, mpiPrint, mpiBarrier, loop resume

time: Decom    MasterSlave
1: 145.04s      146.42s
2:  71.87 s     148.59s
4:  50.3435 s     54.435s
8:  27.573s       23.21s
16: 16.7267s      15.437s
32: 9.356         17.01s
