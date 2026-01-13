# Montecarlo
Montecarlo Integration  
This project implements a Monte Carlo integration algorithm in C++ to calculate definite integrals over N-dimensional hyperspheres and hyperrectangles. The software can use the **muParserX** library to parse mathematical functions at runtime.

## ⚙️ Run Modes

The program offers six different run modes, each designed for specific use cases:

### Mode 1: Parser-Based Integration (Function from File)
Uses the mathematical expression from `function.txt` via the **muParserX** library.

**Use case:** When you want to test custom functions without recompiling.  
**Performance:** Slower due to runtime parsing overhead.  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** The parser reads your function expression and performs Monte Carlo integration using uniform sampling over various geometric domains (hypersphere, hyperrectangle, hypercylinder, etc.).

### Mode 2: Hardcoded Integration (Uniform Distribution)
Uses pre-compiled C++ lambda functions with standard uniform Monte Carlo sampling.

**Use case:** When performance is critical and you're working with fixed functions.  
**Performance:** Fast (no parsing overhead).  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** Executes integration benchmarks across multiple domains using hardcoded test functions with the classic Monte Carlo estimator.

### Mode 3: Hardcoded Integration (Metropolis-Hastings)
Uses pre-compiled C++ lambda functions with Metropolis-Hastings MCMC sampling.

**Use case:** For complex or high-dimensional domains where uniform sampling is inefficient.  
**Performance:** Fast, with better convergence for difficult geometries.  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** Combines volume estimation (Hit-or-Miss) with MH sampling to explore the domain more efficiently, especially useful when the domain has complex boundaries.

### Mode 4: Polytope Integration (Convex Hull)
Integrates over arbitrary convex polytopes defined by user-provided points.

**Use case:** When your integration domain is a convex polytope (e.g., polyhedron in 3D).  
**Performance:** Depends on polytope complexity.  
**Gnuplot:** Not applicable.

**How it works:** Reads vertex coordinates from `points.txt` and facet normals/offsets from `hull.txt` (generated via Qhull). Performs Monte Carlo integration over the resulting convex polytope using the half-space representation.

### Mode 5: Particle Swarm Optimization (PSO)
Runs optimization benchmarks using the Particle Swarm Optimization algorithm.

**Use case:** Finding global minima/maxima of continuous functions.  
**Performance:** Good for smooth, continuous optimization landscapes.  
**Gnuplot:** May generate convergence plots depending on implementation.

**How it works:** Initializes a swarm of particles that explore the search space, updating positions based on their own best solutions and the global best. Tests on standard benchmark functions (Rastrigin, Rosenbrock, etc.).

### Mode 6: Genetic Algorithm (GA)
Runs optimization benchmarks using a Genetic Algorithm.

**Use case:** Finding global minima/maxima, especially for non-smooth or discrete problems.  
**Performance:** Robust for complex fitness landscapes.  
**Gnuplot:** May generate convergence plots depending on implementation.

**How it works:** Evolves a population of candidate solutions through selection, crossover, and mutation operators. Tests on standard benchmark functions to measure convergence and solution quality.

---

### General Notes
- **Seed Control:** You can specify a custom random seed by running: `./montecarlo_1 <seed>`
- **Gnuplot Visualization:** Available for integration modes 1-3 in dimensions 1D, 2D, and 3D only.
- **Closing Plots:** To close all gnuplot windows, run:
  ```bash
  pkill -f gnuplot
  ```

## 📝 How to Write the Function (muParserX)

Before compiling the program, enter the function to integrate in the function.txt file located in the root of the repository.  
Use standard mathematical notation such as:

* **Variables:** Use the variable names defined in your code (e.g., `x`, `y`, `z` or `x[0]`, `x[1]`).
* **Operators:** `+`, `-`, `*`, `/`, `^` (power).
* **Common Functions:** `sin()`, `cos()`, `tan()`, `exp()`, `log()`, `sqrt()`, `abs()`.
* **Constants:** `_pi`, `_e`.

**Example Input:**
```text
sin(x) * exp(y) + (z^2 / 2)
```
## How to set gnuplot on mac

1. Open XQuartx and enter in its terminal:
  ```bash
  xhost +localhost
  ```
2. In container terminal enter:
  ```bash 
  export DISPLAY=host.docker.internal:0
  gnuplot -persist -e "set term x11"
  ```

## 🚀 Build and Run
Ensure you have CMake, a C++ compiler, and muParserX installed.  
From the root of the project:  
  
1. Create the build directory and enter it:
  ```bash
  mkdir build
  ```
2. Generate the build files:
  ```bash
  cmake ..
  ```
3. Compile the project:
  ```bash
  make
  ```
4. Run the executable:
  ```bash
  ./montecarlo_1
  ```

# Mathematical Background

This library implements several Monte Carlo methods for numerical integration and stochastic optimization. Below is the mathematical formulation for the core algorithms used.

## 1. Numerical Integration

### 1.1 Classic Monte Carlo Integration

The simplest Monte Carlo estimator approximates the integral of a function $f$ over a domain $\Omega$ by sampling points uniformly within a bounding box $B$ (where $\Omega \subseteq B$) and evaluating $f$.

$$I = \int_{\Omega} f(\mathbf{x}) \, d\mathbf{x} \approx V_B \cdot \frac{1}{N} \sum_{i=1}^{N} f(\mathbf{x}_i) \cdot \mathbb{1}_{\Omega}(\mathbf{x}_i)$$

Where:

* $V_B$ is the volume of the hyper-rectangle bounding the domain.
* $\mathbf{x}_i$ are $N$ samples drawn uniformly from $B$.
* $\mathbb{1}_{\Omega}(\mathbf{x})$ is the indicator function ($1$ if $\mathbf{x} \in \Omega$, $0$ otherwise).

**Implementation Details:**
The method `integrate` generates random points in the bounding box. If a point falls inside the domain (checked via `domain.isInside(p)`), its contribution is added. Non-domain points effectively contribute 0.

### 1.2 Importance Sampling

Importance Sampling reduces variance by sampling points from a proposal distribution $q(\mathbf{x})$ that ideally resembles the shape of $f(\mathbf{x})$, rather than sampling uniformly.

$$I = \int_{\Omega} f(\mathbf{x}) \, d\mathbf{x} \approx \frac{1}{N} \sum_{i=1}^{N} \frac{f(\mathbf{x}_i)}{q(\mathbf{x}_i)}$$

Where:

* $\mathbf{x}_i$ are samples drawn from the probability density function (PDF) $q(\mathbf{x})$.
* The term $\frac{f(\mathbf{x}_i)}{q(\mathbf{x}_i)}$ is the importance weight.

**Implementation Details:**
The `ISMeanEstimator` computes the mean of the weighted samples. The `integrate_importance` method uses this estimator. Note that for the integration to be correct over the domain volume, the weights must be properly scaled relative to the domain measure.

### 1.3 Metropolis-Hastings Integration

This method is effective for complex domains where simple uniform sampling is inefficient. It separates the problem into two parts: estimating the domain volume and estimating the function mean.

1. **Volume Estimation ($V_{\Omega}$):** Computed using Hit-or-Miss Monte Carlo (see section 3).
2. **Mean Estimation ($\bar{f}$):** Uses the Metropolis-Hastings (MH) algorithm to generate a Markov Chain of samples distributed according to a target density $\pi(\mathbf{x})$ (typically uniform over $\Omega$).

$$I = V_{\Omega} \cdot \bar{f} \quad \text{where} \quad \bar{f} = \frac{1}{N} \sum_{i=1}^{N} f(\mathbf{x}_i)$$

Where $\mathbf{x}_i$ are samples generated by the MH sampler. The sampler accepts a candidate $\mathbf{x}'$ from current state $\mathbf{x}$ with probability:

$$\alpha = \min\left(1, \frac{\pi(\mathbf{x}')}{\pi(\mathbf{x})}\right)$$

(Note: The proposal distribution in the random walk is symmetric, so the Hastings ratio $\frac{q(\mathbf{x}|\mathbf{x}')}{q(\mathbf{x}'|\mathbf{x})}$ cancels out).

**Implementation Details:**
The method `integrate_with_mh` first calls `VolumeEstimatorMC` to find $V_{\Omega}$. Then it runs a `MetropolisHastingsSampler` to gather samples $\mathbf{x}_i$ inside $\Omega$ and computes their average $\bar{f}$.

## 2. Optimization Algorithms

### 2.1 Particle Swarm Optimization (PSO)

PSO optimizes a function by simulating a swarm of particles moving through the search space. Each particle $i$ has a position $\mathbf{x}_i$ and velocity $\mathbf{v}_i$. They are updated based on:

1. **Inertia:** The particle's previous velocity.
2. **Cognitive Component:** The distance to the particle's own best known position ($\mathbf{p}_i$).
3. **Social Component:** The distance to the swarm's global best known position ($\mathbf{g}$).

**Update Equations:**

$$\mathbf{v}_i(t+1) = w \cdot \mathbf{v}_i(t) + c_1 \cdot r_1 \cdot (\mathbf{p}_i - \mathbf{x}_i(t)) + c_2 \cdot r_2 \cdot (\mathbf{g} - \mathbf{x}_i(t))$$

$$\mathbf{x}_i(t+1) = \mathbf{x}_i(t) + \mathbf{v}_i(t+1)$$

Where:

* $w$: Inertia weight.
* $c_1, c_2$: Cognitive and social coefficients.
* $r_1, r_2$: Random numbers in $[0,1]$.

**Implementation Details:**
Found in `PSO.cpp`, the `step()` function applies these updates and enforces boundary constraints (damping velocity if a boundary is hit).

### 2.2 Genetic Algorithm (GA)

GA evolves a population of candidate solutions using biologically inspired operators.

1. **Selection (Tournament):** $k$ individuals are chosen at random, and the best one is selected for reproduction.
2. **Crossover (Uniform):** Two parents swap genes (coordinates) with probability $p_c$ per dimension to create offspring.
3. **Mutation (Gaussian):** Each gene is perturbed by a Gaussian noise $\mathcal{N}(0, \sigma^2)$ with a small probability.
4. **Elitism:** The top $n_e$ best individuals are carried over unchanged to the next generation.

**Implementation Details:**
Found in `GA.cpp`, implementing tournament selection, uniform crossover, and Gaussian mutation within the `step()` loop.

## 3. Volume Estimation

The volume of an integration domain $\Omega$ is estimated using a Hit-or-Miss approach. By enclosing $\Omega$ in a bounding box $B$ of known volume $V_B$:

$$V_{\Omega} \approx V_B \cdot \frac{N_{\text{hits}}}{N_{\text{total}}}$$

The standard error of this estimate is derived from the variance of the Bernoulli distribution (inside/outside):

$$\sigma_V = V_B \cdot \sqrt{\frac{p(1-p)}{N}}$$

Where $p = \frac{N_{\text{hits}}}{N_{\text{total}}}$.

## 📐 Using Qhull for Polytopes

The program also supports Monte Carlo integration over convex polytopes in any dimension.
To define a polytope, you provide a set of points and let Qhull compute its convex hull, including facet normals and offsets.

**1️⃣ Create the input file points.txt**

  The format must be:
  
  ```text
    <num_points> <dim>
    x₁ y₁ z₁
    x₂ y₂ z₂
    ...
  ```
  Example for 3D:
  
  ```text
    10 3
    0.1 0.3 0.5
    0.2 0.8 0.4
    ...
  ```

**2️⃣ Compute the convex hull (normals + offsets)**

  Use Qhull with the n option to output facet normals and plane offsets:
  
  ```bash
    module load qhull
    qhull Qt n < points.txt > hull.txt
  ```
  Where:
  	•	Qt → triangulates the hull
  	•	n  → prints one facet normal per line, followed by its offset d
  	•	Qhull outputs hyperplanes in the form
  
  n · x + d = 0
  
  and the program internally converts them to
  
  n' · x ≤ b
  
  which defines the half-spaces forming the convex polytope.

**3️⃣ Run the program in polytope mode**
