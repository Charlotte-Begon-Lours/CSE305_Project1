# CSE305_Project1
## CSE305 Project: N-body simulations
Elisabeth Le Leslé - Fanny Scherer - Charlotte Bégon-Lours

## Introduction 
Here's a modified version of your LaTeX content rewritten in a clear and professional **README.md** format, suitable for a GitHub or project repository:

---

# N-Body Simulation in 2D

## Overview

This project implements and analyzes an **N-body simulation system** for modeling interactions between bodies in 2D space. This is a fundamental problem in physics, involving the simulation of multiple objects—such as planets, satellites, and stars—interacting through gravitational forces.

The goal of this project is to:

* Implement multiple simulation methods.
* Explore parallelization strategies.
* Provide results for analyzing the **efficiency**, **accuracy**, and **scalability** of different algorithmic approaches across varying thread configurations and problem sizes.

Given **N bodies** with initial positions, velocities, and masses, the objective is to simulate their **gravitational interactions over time** using **Newton's law of motion**. According to the **universal law of gravitation**, each body experiences forces from all other bodies, with force inversely proportional to the square of the distance between them.

Each body `bᵢ` is defined by:

* `mᵢ`: Mass
* `(xᵢ, yᵢ)`: Position in 2D space
* `(vxᵢ, vyᵢ)`: Velocity in 2D space

Two bodies `bᵢ` and `bⱼ` exert a mutual force calculated as:

```
fᵢⱼ = G * (mᵢ * mⱼ) / ((xᵢ - xⱼ)² + (yᵢ - yⱼ)²)
```

Where:

* `G = 6.67430 × 10⁻¹¹ m³·kg⁻¹·s⁻²` is the gravitational constant.



## Methods implemented

- Sequential simulation of the N-body system
- Parallelization of the update step
- Parallelization of the force computations using mutex
- Parallelization of the force computations without mutex
- Barnes-Hut
- FMM


## How to run the code

### When compiling:

$ make ./simple_approach_test

### When running: 

$ ./simple_approach_test

### Run the Python file in venv environment to generate the pngs for the simulations:

(venv) $ python make_frames.py

### Back in terminal, an example of generating a gif from the png frames using Magick++:

$ magick convert -delay 20 -loop 0 frames1/frame_*.png simulation1.gif

### To run the plots of the different simulation runtimes w.r.t. the different threads:

First modify the main() in the simple_approach_test.cpp.

$ make clean 

$ make ./simple_approach_test

(venv) $ python benchmark_pipeline.py


### To run the plots of the different simulation accuracies w.r.t the number of bodies:

Modify the main() in the simple_approach_test.cpp.

$ make clean 

$ make ./simple_approach_test

(venv) $ python accuracy_pipeline.py





