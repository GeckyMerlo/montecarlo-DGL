# Project Roadmap & TODOs

## Next Imminent Steps
- [ ] **Benchmarking Suite:** Create a script containing different premade domains and shapes (hyperrectangles and hyperspheres) to test against various mathematical functions automatically.
- [ ] **Enhanced OpenMP:** Improve parallelization strategies, may also extending it to handle the concurrent computation of multiple integrals on multiple domains.
- [ ] **Parser Performance Analysis:** - Quantify the performance gap between **Hardcoded functions** (compile-time) vs. **Parsed functions** (runtime).
    - Compare **muParser** (simpler, faster) vs. **muParserX** (more complex features) to determine the trade-off between speed and functionality.

## More Advanced Implementations

### 1. Advanced Integration Algorithms
- [ ] **Explore Variance Reduction:** Investigate techniques beyond simple Monte Carlo (e.g., Stratified Sampling, Importance Sampling).
- [ ] **Metropolis-Hastings:** Implement the Metropolis algorithm or other MCMC methods to handle complex distributions more efficiently.
- [ ] **Quasi-Monte Carlo:** Experiment with low-discrepancy sequences (e.g., Sobol sequences) to compare convergence rates against pseudo-random numbers.

### 2. Geometric Domain Expansion
- [ ] **Polygons & Polyhedra:** Extend integration support from simple hyperspheres/boxes to arbitrary 2D polygons and 3D polyhedra.
- [ ] **N-Dimensional Polytopes:** Generalize the geometry to support N-dimensional polytopes.
- [ ] **Convexity Handling:**
    - [ ] Implement specific optimizations for **convex** shapes.
    - [ ] Develop rejection sampling or triangulation strategies for **non-convex** domains. **(HARD)**

### 3. Parallelization
- [ ] **Multi-threading (CPU):** Implement **OpenMP** directives to parallelize the sample generation loop across available CPU cores. **(WIP)**
- [ ] **GPU Acceleration (CUDA):** Port the core Monte Carlo kernel to **CUDA** to leverage massive parallelism on NVIDIA GPUs.
- [ ] **Benchmarking:** Create a comparison report of execution times: Serial vs. OpenMP vs. CUDA.

### 4. Practical Applications & Case Studies
- [ ] **Literature Review:** Analyze `Quinn_CH10.pdf` (located in the repo) to identify standard parallel computing problems solvable via Monte Carlo.
- [ ] **Physics Implementation:** Simulate physical phenomena (e.g., heat transfer, particle diffusion) based on the text.
- [ ] **Mathematical Modeling:** Apply the integrator to higher-dimensional mathematical problems suggested in the reference material.
