# Montecarlo
Montecarlo Integration  
This project implements a Monte Carlo integration algorithm in C++ to calculate definite integrals over N-dimensional hyperspheres and hyperrectangles. The software can use the **muParserX** library to parse mathematical functions at runtime.

## ⚙️ Configuration
When the program runs, it will ask you if you want to use the function written in function.txt or use the hardcoded functions inside.  
After that it will ask if you want to plot (Remember: it only works for the following dimensions: 1D,2D,3D)  
To close all the windows use the following command from terminal:
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

# Montecarlo
Montecarlo Integration  
This project implements a Monte Carlo integration algorithm in C++ to calculate definite integrals over N-dimensional hyperspheres and hyperrectangles. The software can use the **muParserX** library to parse mathematical functions at runtime.

## ⚙️ Configuration
When the program runs, it will ask you if you want to use the function written in function.txt or use the hardcoded functions inside.  
After that it will ask if you want to plot (Remember: it only works for the following dimensions: 1D,2D,3D)  
To close all the windows use the following command from terminal:
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
