# BS/GBS Using Tensor Networks

## Project Overview

This is a small part of a project conducted for 'Computation, Learning, and Physics, 2025 Spring'. The code is a crude implementation based on the two following papers.

- Cilluffo, D., Lorenzoni, N., & Plenio, M. B. (2023). Simulating Gaussian Boson Sampling with Tensor Networks in the Heisenberg picture. In arXiv [quant-ph]. arXiv. http://arxiv.org/abs/2305.11215
- Wu, Y.-H., Wang, L., & Tu, H.-H. (2020). Tensor network representations of parton wave functions. Physical Review Letters, 124(24), 246401.

The scalability of this code is not great, and moreover in small dimensions it's more beneficial to calculate the probability using the Hafnian (see resources such as StrawberryFields), so this implementation is more or less pepedagogical. In future works, GPU acceleration may be able to help.

## Description

- The codes inside `./tensor` and `./utils`, and part of `startup.m` are written by S.-S. B. Lee., and is contained in https://github.com/ssblee/TensorNetworks2022.
- The functions for calculating BS/GBS probabilities are located in `./linearoptics`.
- The DMFT-based method for calculating the MPO-MPS contraction is based on 1-site DMFT. This will get extended to 2-site DMFT in future works.
