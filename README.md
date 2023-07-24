# HPEZ

## Introduction

This is the source code of the HPEZ compressor introduced in the paper: High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation.

## Dependencies

Please Install the following dependencies before running the artiact evaluation experiments:

* cmake>=3.13
* gcc>=6.0

## 3rd party libraries/tools

* Zstd >= 1.3.5 (https://facebook.github.io/zstd/). Not mandatory to be mannually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include. A Cmake version >= 3.13.0 is needed. 
Before you proceed to the following evaluations, please add the installation path of HPEZ to your system path so that you can directly run qoz command in your machine for further evaluations.

## Single compression/decompression testing Examples

You can use the executable 'hpez' command to do the compression/decompression (the input data should be float or double binary files). Just run "hpez" command without any argument to check the instructions for its arguments.
For the convenience of tests, the hpez executable includes the SZ3.1 compression, QoZ 1.1 compression, and 3 optimization levels of HPEZ compression. In the command:
* Not containing -q argument or -q 0: SZ3.1 compression.
* Containing -q 1: QoZ 1.1 compression.
* Containing -q 4: Full HPEZ compression (for the results reported in the paper).
* Containing -q 2 or -q 3: 2 intermediate optimization levels of HPEZ compression (having faster speeds but slightly worse rate-distortion).

Notice: the intergrated SZ3.1 and QoZ 1.1 in HPEZ has already leveraged the Fast-varying-first interpolation (proposed in our paper), therefore their compression ratios are sometimes higher than the original public released versions of SZ3.1 and QoZ 1.1.

## Dataset Source 

Please download test datasets from: https://sdrbench.github.io/. 

