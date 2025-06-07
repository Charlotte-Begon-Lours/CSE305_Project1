# CSE305_Project1
## CSE305 Project: N-body simulations
Elisabeth Le Leslé - Fanny Scherer - Charlotte Bégon-Lours

## Introduction 
The purpose of this project was to implement and analyse a N-body simulation system for modeling interactions between bodies in 2D space. This challenge takes major importance in physics, involving the simulation of multiple objects interacting through gravitational forces such as planets, satellites and stars. \\
\\
Hence, the goal of this project is to implements multiple methods and parallelization strategies for the simulation, providing valuable results for analyzing the efficiency of those methods as well and the accuracy and scalability of different algorithmic approaches across various thread configurations and problem sizes. \\

Given N bodies with initial positions, velocities, and mass, we want to simulate their gravitational interactions over time using Newton's law of motion. From the universal law of gravitation, each body experiences forces from all other bodies inversely proportional to the distance that separates them.\\
\\
Each body $b_{i}$ is characterized by:
\begin{itemize}
    \item a mass $m_{i}$
    \item a position in 2D space $(x_{i}, y_{i})$
    \item a velocity in 2D space $(vx_{i}, vy_{i})$
\end{itemize}
Two bodies $b_{i}$ and $b_{j}$ act on each other through a force given by:
\begin{align}
f_{i,j} = G \times \frac{m_i m_j}{(x_i - x_j)^2 + (y_i - y_j)^2}
\end{align}
where \( G = 6.67430 \times 10^{-11} \, m^{3} \cdot kg^{-1} \cdot s^{-2} \) is the gravitational constant.
\\


## Methods implemented

- Sequential simulation of the N-body system
- Parallelization of the update step
- Parallelization of the force computations using mutex
- Parallelization of the force computations without mutex
- Barnes-Hut
- FMM


## How to run the code

When compiling: $ make ./simple_approach_test
When running: $ ./simple_approach_test

Run the Python file in venv environment to generate the pngs for the simulations:
(venv) $ python make_frames.py

Back in terminal, an example of generating a gif from the png frames using Magick++:
$ magick convert -delay 20 -loop 0 frames1/frame_*.png simulation1.gif

To run the plots of the different simulation runtimes w.r.t. the different threads:
First modify the main() in the simple_approach_test.cpp.
$ make clean 
$ make ./simple_approach_test
(venv) $ python benchmark_pipeline.py

To run the plots of the different simulation accuracies w.r.t the number of bodies:
Modify the main() in the simple_approach_test.cpp.
$ make clean 
$ make ./simple_approach_test
(venv) $ python accuracy_pipeline.py





