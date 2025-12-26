# TRFx: Unlocking Genome-Wide Tandem Repeat Analysis via a Hybrid CPU-GPU Acceleration of Tandem Repeat Finder
50x faster TRF for modern genomics | Strict output consistency | GPU acceleration | Multi-threading

TRFx is a highly optimized, parallel version of the legendary Tandem Repeat Finder ([TRF4.10](https://github.com/liyanhui607/TRF4.10)), designed to overcome the computational bottleneck of the original algorithm in the era of third-generation sequencing (TGS) data. It achieves 50x speedups​ on standard servers, reducing the analysis of human-scale TGS datasets from weeks down to hours, while maintaining strict output consistency with the original TRF.


## Environment:
Operating System: The server ran on the Ubuntu 20.04.6 LTS operating system.

Compiler: All code was compiled using GCC version 9.4.0.

GPU Computing: GPU-accelerated computations were supported by the CUDA toolkit version 12.8.



## One-Command Demo (Recommended for First-Time Users)
```
git clone https://github.com/liyanhui607/TRFx.git

cd TRFx/test

bash cmd 

```


This will:
✅ Automatically compile TRFx
✅ Run analysis on example data (FASTA)
✅ Use 20 CPU threads
✅ Show detailed performance statistics
✅ Save results to output.txt


## Installation
```
git clone https://github.com/your_username/TRFx.git

cd TRFx

make

```

## Usage examples
Use all default parameters except threads

`./trfx input.fasta -t 20`

 Quick analysis with minimal parameters
 
`./trfx reads.fastq -t 8 > results.txt`

Custom parameters for specific analysis

`./trfx genome.fa -a 2 -b 7 -d 7 -m 80 -i 10 -s 50 -p 2000 -t 20 > output.txt`



## Parameters

Default: 'trfx inputFile -a 2 -b 7 -d 7 -m 80 -i 10 -s 50 -p 2000 -t 3'

Default is good in most time , So simply use: ./trfx inputFile -t 8 (number of threads)

Where: (all weights, penalties, and scores are positive)

  File = sequences input file
  
  Match  = matching weight [2]
  
  Mismatch  = mismatching penalty [7]
  
  Delta = indel penalty [7]
  
  PM = match probability (whole number) [80]
  
  PI = indel probability (whole number) [10]
  
  Minscore = minimum alignment score to report [50]
  
  MaxPeriod = maximum period size to report [2000]
  
  [options] = one or more of the following:
  
  -t INT     number of threads [3]
  
  -V         show version number


## Frequently Asked Questions
Q: Do I need to adjust all parameters?

A: No. Default parameters work well for most applications. Just specify your input file and number of threads.

Q: How many threads should I use?

A: Typically use the number of CPU cores available. For modern servers, 16-32 threads is common.

Q: Is GPU required?

A: No. Significant speedups (20-30x) are achievable with CPU threads alone. GPU provides additional boost for long sequences.

Q: When should I adjust parameters?

A: Only if you need specific sensitivity tuning for specialized analyses.

Q: Where could I find the original TRF4.10?

A: [TRF4.10](https://github.com/liyanhui607/TRF4.10) can be found here: https://github.com/liyanhui607/TRF4.10.



## Technical Innovations
TRFx accelerates TRF through three key optimizations:

1. Multi-threaded Pipeline Architecture
   
I/O-Compute Overlap: Parallel data reading, processing, and writing

Work Stealing: Dynamic load balancing across threads

Bulk Data Loading: Reads large blocks (10GB+) to minimize disk I/O


2. CPU Optimizations
   
Memory Access Patterns: 16 separate arrays for 2-mer processing

Dynamic Memory Reshaping: On-demand S-array column allocation

Modulo Operation Replacement: Conditional logic for expensive divisions


3. Hybrid CPU/GPU Acceleration
Smart Workload Distribution: 

Sequences >4000bp automatically routed to GPU

Shared Memory Optimization: Local histogram accumulation reduces global memory contention

Adaptive Block Sizing: CUDA occupancy API for optimal GPU utilization


## Citation
If you use TRFx in your research, please cite:

Yan-Hui Li, Li Fang, Yuan Zhou. TRFx: Accelerating Tandem Repeat Finder with Multithreading and GPU for Third-Generation Sequencing Data. [Journal Name, Volume, Pages, Year].



## Contributing
We welcome contributions! Please feel free to submit issues, feature requests, or pull requests.


## License
TRFx is an optimized, parallel derivative of the Tandem Repeats Finder (TRF) by Benson et al. In compliance with the original work's license, TRFx is released under the GNU Affero General Public License v3.0 (AGPL-3.0).

This license permits use, modification, and distribution for any purpose, including commercial use, provided all distribution conditions of the AGPL-3.0 are met.

For inquiries regarding alternative licensing arrangements (e.g., for proprietary integration), please contact the corresponding author.

The GitHub repository contains the full source code, build instructions, documentation, and examples to ensure reproducibility and foster community adoption.
