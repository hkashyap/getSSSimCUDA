# getSSSimCUDA
CUDA C implementation of getSSSim gene coexpression similarity measure

# Tested using:
Platform: Ubuntu 16.04, gcc 5.4.0, CUDA 7.5

Graphics card: NVIDIA GeForce GTX 980

# To use:
1. The provided Makefile compiles the allGenes.cu, which is the main CUDA implmentation of getSSSim
2. For users without a GPU, we also provide a CPU only implementation of getSSSim. You can compile it using g++ (without CUDA API).
3. Additionally, we provide the version of the getSSSim used in our 2016 Nature Scientific Reports paper.

If you use this code in your research, please cite our paper:

Kakati, Tulika, Hirak Kashyap, and Dhruba K. Bhattacharyya. "THD-module extractor: an application for CEN module extraction and interesting gene identification for Alzheimer’s disease." Scientific reports 6 (2016): 38046.
