#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <string>
#include <sstream>

const double G = 6.67430e-11; // Gravity constant

// Structure to represent a body in 2D space: we use a struct instead of a class as we want to make it publicly accessible in the whole code
struct Body {
    double mass;     // of the body
    double x, y;     // position coordinates
    double vx, vy;   // velocity components
    double fx, fy;   // Force components

}
