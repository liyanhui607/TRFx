# TRFx: Unlocking Genome-Wide Tandem Repeat Analysis on Third-Generation Sequencing Data
**40x faster | Bit-identical output consistency | GPU acceleration | Multi-threading**

TRFx is a highly optimized, parallel version of the legendary Tandem Repeat Finder ([TRF4.10](https://github.com/Benson-Genomics-Lab/TRF)), designed to overcome the computational bottleneck of the original algorithm with third-generation sequencing (TGS) data. It achieves more than 40x speedups​ on standard servers (20 cores, 20 threads), reducing the analysis of human-scale TGS datasets from weeks down to hours, while maintaining bit-identical​​ output with the original TRF and ([TRF-mod](https://github.com/lh3/TRF-mod)).




## Key Features
**Unprecedented Speed**: 40x overall acceleration achieved through synergistic optimization.

**Bit-identical Output Consistency**: All output formats (dat, mask, ngs) are bit-identical​ with the original TRF. The default BED output is bit-identical with TRF-mod.

**Hybrid CPU/GPU Architecture**: Intelligent workload distribution automatically routes long sequences  to the GPU for additional acceleration (typically 1.25x-1.4x improvement).

**Minimal Parameter Tuning**: Default parameters work optimally for most applications.




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

This command will:

1. Navigate to the project root directory

2. Clean previous compilation files (make clean)

3. Compile TRFx (make)

4. Return to the testdirectory

5. Run TRFx on example data input.fasta with parameters -a 2 -b 7 -g 7 -k 80 -i 10 -s 30 -p 2000 -t 8

6. Redirect BED format result to output.bed




## Installation
```
git clone https://github.com/liyanhui607/TRFx.git

cd TRFx

make

```



## Usage
**Basic Syntax**

`trfx [options] <in.fa>`


**Examples**

```
# Basic usage with default parameters (8 threads)

./trfx your_sequence.fasta

# Specify thread count for better performance

./trfx your_sequence.fasta -t 20

# Comprehensive analysis with all output formats

./trfx -a 2 -b 7 -g 7 -k 80 -i 10 -s 30 -p 2000 -t 8 -m -d -n your_sequence.fasta > your_sequence.fasta.2.7.7.80.10.30.2000.bed

```

Simply replace your_sequence.fasta with your actual FASTA file path.​ The comprehensive command generates all four output formats simultaneously.




## Parameters
**Core Parameters**

| Parameter | Description | Default Value |
|:---|:---|:---|
| `-a INT` | Match weight (matching score) | 2 |
| `-b INT` | Mismatch penalty | 7 |
| `-g INT` | Indel penalty | 7 |
| `-k INT` | Match probability (75 or 80) | 80 |
| `-i INT` | Indel probability (10 or 20) | 10 |
| `-s INT` | Minimum alignment score to report (30 or 50)  | 30 |
| `-p INT` | Maximum period size to report [1-2000] | 500 |
| `-l INT` | Maximum TR length expected (in millions) | 12 |
| `-t INT` | Number of threads | 8 |


**Output Format Options**

-n: Output in TRF NGS format.

-d: Output in TRF dat format.

-m: Output masked sequence file.

-r: Don't eliminate redundancy.

-v: Print version information.




##  Output Files
TRFx generates multiple output formats to support different downstream analyses. Output files follow this naming pattern:

`<input_filename>.<a>.<b>.<g>.<k>.<i>.<s>.<p>.<extension>\n`


### Output File Specifications

**1. BED Format (.bed):**

Default output to stdout (typically redirected to a file).

Fully compatible with TRF-mod, using standard BED format suitable for genome browsers and various genomic analysis tools.

Format: ctg start end period copyNum fracMatch fracGap score entropy pattern

Note: Pattern length may differ from period in BED format.

**2. TRF Dat Format (.dat):**

Generated with the -d option.

Bit-identical​ to the original TRF's dat output, containing detailed alignment information. 

**3. TRF Mask Format (.mask):**

Generated with the -m option.

Bit-identical​ to the original TRF's mask output, providing masked (tagged) sequences. 

**4. TRF NGS Format (.ngs):**

Generated with the -n option.

Bit-identical​ to the original TRF's ngs output, offering comprehensive statistics for next-generation sequencing analysis. 




## Technical Architecture
TRFx accelerates TRF through three key optimizations:


**1. Multi-threaded Pipeline Architecture:** Implements I/O-compute overlap with dynamic load balancing and bulk data loading to minimize disk I/O.

**2. CPU Optimizations:** Features enhanced memory access patterns, dynamic memory reshaping, and replacement of expensive modulo operations with conditional logic.

**3. Hybrid CPU/GPU Acceleration:** Employs smart workload distribution with shared memory optimization and adaptive block sizing using CUDA occupancy API for optimal GPU utilization.




##  Frequently Asked Questions
**Q: How is the 40× speedup calculated ?**

**A**:​ The 40× acceleration is measured on a 20-core server running 20 threads compared to single-threaded execution of the optimized TRFx code. It represents the combined effect of single-threaded optimizations (2.55-3.04×) and excellent parallel scaling (16.45-16.76× on 20 threads).

**Q: Which output formats should I use ?**

**A**:​ Use the default BED format for compatibility with genomic tools and browsers (like TRF-mod). Use the -d, -m, and -n formats for backward compatibility with existing pipelines that rely on TRF's native outputs.

**Q: Is GPU required for running TRFx ?**

**A**:​ No. Significant speedups (20-30×) are achievable with CPU threads alone. GPU provides an additional boost (typically 1.25-1.4×).

**Q: When should I adjust parameters from defaults?**

**A**:​ Only for specialized analyses requiring specific sensitivity tuning. Default parameters are optimized for most applications.

**Q: Why is "TRF" capitalized and the "x" in lowercase in TRFx?**

**​A​**:​ The capitalization signifies respect for the original Tandem Repeat Finder (TRF) software upon which TRFx is built. 




## Citation
If you use TRFx in your research, please cite:

Yan-Hui Li, Li Fang, Yuan Zhou. TRFx: Unlocking Genome-Wide Tandem Repeat Analysis on Third-Generation Sequencing Data. [Journal Name, Volume, Pages, Year].




## Contributing
We welcome contributions! Please feel free to submit issues, feature requests, or pull requests.




## License
TRFx is an optimized, parallel derivative of the Tandem Repeats Finder (TRF) by Benson et al. In compliance with the original work's license, TRFx is released under the GNU Affero General Public License v3.0 (AGPL-3.0).

This license permits use, modification, and distribution for any purpose, including commercial use, provided all distribution conditions of the AGPL-3.0 are met.

For inquiries regarding alternative licensing arrangements (e.g., for proprietary integration), please contact the corresponding author.


