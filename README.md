# getSSSimCUDA
CUDA C implementation of getSSSim gene coexpression similarity measure

# Tested using:
Platform: Ubuntu 14.04, gcc 4.8.4, CUDA 7.5

Graphics card: NVIDIA GeForce GTX 980

# To use:
1. The provided Makefile compiles the allGenes.cu, which is the main CUDA implmentation of getSSSim
2. For users without a GPU, we also provide a non cuda implementation of getSSSim, which is still compiles using nvcc, as it uses some functions from the CUDA library for storekeeping. Can be converted to use standard C library functions to compile uisng GCC alone (without CUDA API).
3. We also provide the version of the getSSSim used in our 2016 Nature Scientific Reports paper.

If you use this library in your research, please cite our paper:

Kakati, Tulika, Hirak Kashyap, and Dhruba K. Bhattacharyya. "THD-module extractor: an application for CEN module extraction and interesting gene identification for Alzheimer’s disease." Scientific reports 6 (2016): 38046.
