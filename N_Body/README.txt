Tomasz Michalik MPCS 51087 - nBody problem

**notes: had issues using all files together. Submitting full files but MPI testing file is : pipeline.cpp or PIPE exe. Copied all methods needed for MPI
Got weird error in full version where NumberIterations went to 0 randomly so wasnt able to output full file. Otherwise all
should be implemented

implementations: (home comp)
MPI
Serial : 13.75s avg
OpenMP 2:  11.74s    4:  10.02s  8: 9.92s  16: 9.65s

Not completely same when running diff - though animations are similar (see SerCorrect vs ParallelCorrect mp4s).
Data is in (ParCorrect & SerCorrect txt files)
Due to randomness i dont believe i should get the same random numbers as serial version. May be due to my method of
having master send pieces to all ranks

*Swirl.mp4 file showing modded galaxy collision for simulation

midway:
hybrid vs MPI test part 1
mpi: 5.29sec

scaling (midway)
1 node   2 node    4 node       8 node
179 sec  92.2 sec   50.25sec    27.08

perf test(vesta):
with printing 32 nodes MPI full : 29.64s
without printing routine: : 29.18s (suprisingly low when cutting out all printing)

scaling( 512,800 bodies)
64:   375.3 sec
128:  206.5 sec
256:  121.3 sec
512:  93.04 sec
1024: 88.95 sec

****Hybrid implementation - kept gettign a segfault in random times, the Iteration count which i use to create my global size array fails and show 0 ... was not able to resolve

N Bodies =                     102400
Timestep dt =                  2.000e-01
Number of Timesteps =          102400
Number of Threads per Rank =   16
***Outer loop iter:0 counter=102400 rank0
Global: 1024000, rank2,total:1024000 startrow:51200, step:0 local 25600 iter:0
Global: 1024000, rank1,total:1024000 startrow:25600, step:0 local 25600 iter:0
Global: 1024000, rank3,total:1024000 startrow:76800, step:0 local 25600 iter:0
Global: 1024000, rank0,total:1024000 startrow:0, step:0 local 25600 iter:0
***Outer loop iter:1 counter=102400 rank0
Global: 0, rank3,total:0 startrow:76800, step:0 local 25600 iter:1
Global: 0, rank2,total:0 startrow:51200, step:0 local 25600 iter:1
Global: 0, rank1,total:0 startrow:25600, step:0 local 25600 iter:1
Global: 0, rank0,total:0 startrow:0, step:0 local 25600 iter:1