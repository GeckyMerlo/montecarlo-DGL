# Montecarlo Integration
This project implements a Monte Carlo integration algorithm in C++ to calculate definite integrals over N-dimensional hyperspheres and hyperrectangles. The software can use the **muParserX** library to parse mathematical functions at runtime.

## 📚 Documentation

Full API documentation is available via Doxygen:

```bash
# Generate documentation (requires doxygen and graphviz)
doxygen Doxyfile

# Open in browser
open docs/html/index.html  # macOS
xdg-open docs/html/index.html  # Linux
```

See [DOXYGEN_GUIDE.md](DOXYGEN_GUIDE.md) for detailed instructions.

## ⚙️ Run Modes

The project offers several run modes, divided into two main executables (`montecarlo_1` for benchmarks/tests and `drone_optimization` for the specific application):

### Mode 1: Parser-Based Integration (Function from File)
*Executable: `montecarlo_1`* Uses the mathematical expression from `function.txt` via the **muParserX** library.

**Use case:** When you want to test custom functions without recompiling.  
**Performance:** Slower due to runtime parsing overhead.  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** The parser reads your function expression and performs Monte Carlo integration using uniform sampling over various geometric domains (hypersphere, hyperrectangle, hypercylinder, etc.).

### Mode 2: Hardcoded Integration (Uniform Distribution)
*Executable: `montecarlo_1`* Uses pre-compiled C++ lambda functions with standard uniform Monte Carlo sampling.

**Use case:** When performance is critical and you're working with fixed functions.  
**Performance:** Fast (no parsing overhead).  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** Executes integration benchmarks across multiple domains using hardcoded test functions with the classic Monte Carlo estimator.

### Mode 3: Hardcoded Integration (Metropolis-Hastings)
*Executable: `montecarlo_1`* Uses pre-compiled C++ lambda functions with Metropolis-Hastings MCMC sampling.

**Use case:** For complex or high-dimensional domains where uniform sampling is inefficient.  
**Performance:** Fast, with better convergence for difficult geometries.  
**Gnuplot:** Supports visualization (1D, 2D, 3D only).

**How it works:** Combines volume estimation (Hit-or-Miss) with MH sampling to explore the domain more efficiently, especially useful when the domain has complex boundaries.

### Mode 4: Polytope Integration (Convex Hull)
*Executable: `montecarlo_1`* Integrates over arbitrary convex polytopes defined by user-provided points.

**Use case:** When your integration domain is a convex polytope (e.g., polyhedron in 3D).  
**Performance:** Depends on polytope complexity.  
**Gnuplot:** Not applicable.

**How it works:** Reads vertex coordinates from `points.txt` and facet normals/offsets from `hull.txt` (generated via Qhull). Performs Monte Carlo integration over the resulting convex polytope using the half-space representation.

### Mode 5: Particle Swarm Optimization (PSO)
*Executable: `montecarlo_1`* Runs optimization benchmarks using the Particle Swarm Optimization algorithm.

**Use case:** Finding global minima/maxima of continuous functions.  
**Performance:** Good for smooth, continuous optimization landscapes.  
**Gnuplot:** May generate convergence plots depending on implementation.

**How it works:** Initializes a swarm of particles that explore the search space, updating positions based on their own best solutions and the global best. Tests on standard benchmark functions (Rastrigin, Rosenbrock, etc.).

### Mode 6: Genetic Algorithm (GA)
*Executable: `montecarlo_1`* Runs optimization benchmarks using a Genetic Algorithm.

**Use case:** Finding global minima/maxima, especially for non-smooth or discrete problems.  
**Performance:** Robust for complex fitness landscapes.  
**Gnuplot:** May generate convergence plots depending on implementation.

**How it works:** Evolves a population of candidate solutions through selection, crossover, and mutation operators. Tests on standard benchmark functions to measure convergence and solution quality.

### Mode 7: Drone Arm Center of Mass Optimization
*Executable: `drone_optimization`* Specialized application combining Monte Carlo integration with PSO for geometric optimization.

**Use case:** Optimizing hole placement and size in a drone arm to achieve target center of mass.  
**Performance:** Parallelized PSO with OpenMP for efficient multi-core execution.  
**Output:** High-precision verification with 1M samples, auto-generated 3D visualization script.

**How it works:** 
- Defines a complex 3D domain (rectangular arm + cylindrical motor + optional polytope cabin)
- Uses PSO to find optimal hole position [x, y, z] and radius that shifts center of mass to target (1.0, 0, 0)
- Fast optimization with 20,000 samples, followed by high-precision verification with 1,000,000 samples
- Exports sampled geometry to `drone_frames/drone_domain.txt`
- Auto-generates `visualize_drone.gp` script for immediate 3D visualization
- Parallelized particle evaluation with thread-safe RNG management

**Key Findings:**  
The optimization successfully balances the drone arm (Target X: 1.0, Result X: ~0.99) within the physical constraints of the geometry. Crucially, the **Hybrid Verification** (Monte Carlo Baseline - Analytic Hole) matches the pure Monte Carlo re-verification with a delta of **< 0.1%**, confirming the high fidelity and robustness of the stochastic solver.

**Visualization:**  
After running, the script automatically generates `visualize_drone.gp`:
```bash
gnuplot -persist visualize_drone.gp  # Visualize the drone geometry
```

**Images of the domain:**  
*(Note: The images below are from a similar simulation setup to clearly illustrate the hole placement; specific dimensions in the current build may vary slightly.)*

<p float="left">
<img src="drone_assets/drone_domain_view.png" width="48%" alt="Drone domain view" />
<img src="drone_assets/drone_domain_hole.png" width="48%" alt="Drone domain hole detail" />
</p>

**Command line:**
```bash
./drone_optimization [seed] [num_threads]
```
- `seed`: Random seed (optional, default: 12345; use `-` to keep default when specifying threads)
- `num_threads`: Number of threads (optional, default: max available; 0 = sequential for performance comparison)

**Examples:**
```bash
./drone_optimization                 # Run with default seed and max threads
./drone_optimization 42 4            # Run with seed 42 using 4 threads
./drone_optimization - 4             # Run with default seed using 4 threads
./drone_optimization 42 0            # Run with seed 42 sequentially (single thread)
```

### Mode 8: Wind Farm Layout Optimization
*Executable: `wind_farm_simulator`* Optimizes wind turbine placement using hybrid Metropolis-Hastings Monte Carlo integration combined with PSO and Genetic Algorithm.

**Use case:** Finding the optimal layout of wind turbines in a farm to maximize power generation while respecting minimum distance constraints.  
**Performance:** Parallel MH integration with OpenMP thread-safe RNG and optimizers running PSO vs GA comparison.  
**Output:** Optimized turbine positions, convergence plots, and wind farm layout visualizations.

**How it works:**
- Models a 1000m × 1000m wind farm with 15 turbines
- Uses Weibull wind speed distribution and wake effects between turbines
- Employs Metropolis-Hastings integration to estimate average farm power given a turbine layout
- Optimizes turbine positions (x, y for each turbine) using both PSO and GA
- Enforces minimum 50m distance between turbines via penalty function
- Accounts for wind speed reduction downstream of turbines (wake effect)
- Models turbine power output based on wind speed, rotor area, air density, and power coefficient

**Physical Model:**
- **Weibull Wind Distribution:** Wind speeds follow a Weibull distribution (shape k=2.0, scale λ=8.0 m/s)
- **Wake Effect:** Downstream turbines experience reduced wind speed based on distance from upwind turbines
- **Power Model:** Turbine power output: $P = 0.5 \times \rho \times A \times C_p \times v^3$
  - $\rho$ = air density (1.225 kg/m³)
  - $A$ = rotor area (π × 25² m²)
  - $C_p$ = power coefficient (0.4)
  - $v$ = effective wind speed

**Algorithm Configuration:**
- **MH Integration:** 1500 samples, 400 burn-in, thinning factor 2
- **PSO:** 60 particles, 150 iterations, inertia=0.6, cognitive=1.8, social=2.0
- **GA:** 80 population, 200 generations, tournament_k=3, crossover=0.9, mutation=0.1

**Command line:**
```bash
./wind_farm_simulator
```

**Output Files:**
- `results_pso.dat` - Optimized turbine positions from PSO
- `results_ga.dat` - Optimized turbine positions from GA
- `plot_pso.gp` - Gnuplot script for PSO layout visualization
- `plot_ga.gp` - Gnuplot script for GA layout visualization
- `wind_farm_layout_pso.png` - PSO optimization result image
- `wind_farm_layout_ga.png` - GA optimization result image

**Visualization:**
```bash
gnuplot plot_pso.gp  # Visualize PSO optimized wind farm layout
gnuplot plot_ga.gp   # Visualize GA optimized wind farm layout
```

### General Notes
- **Seed Control:** You can specify a custom random seed by running: `./montecarlo_1 <seed>`
- **Parallelization:** PSO and drone_optimization support OpenMP parallelization. Control with:
  - `./montecarlo_1 <seed> <num_threads>` for montecarlo_1
  - `./drone_optimization [seed|-] [num_threads]` for drone optimization
  - Use `0` for sequential execution, `N > 0` for N threads
  - Use `-` as seed placeholder to keep default seed when specifying threads
- **Output Organization:** Files are automatically saved to subdirectories (e.g., `drone_frames/` for the drone geometry outputs).
- **Visualization:**
  - Benchmark animations: `gnuplot run_pso.gp`, `run_pso_3d.gp`, `run_ga.gp`, `run_ga_3d.gp`
  - Drone geometry: `gnuplot -persist visualize_drone.gp` (auto-generated after running drone_optimization)
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
