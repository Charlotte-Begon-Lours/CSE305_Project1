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
const double NOT_ZERO = 1e-9; // constant to avoid division by zero

double sqr(int x) {
    return x*x;
}

// Structure to represent a body in 2D space: we use a struct instead of a class as we want to make it publicly accessible in the whole code
struct Body {
    double mass;     // of the body
    double x, y;     // position coordinates
    double vx, vy;   // velocity components
    double fx, fy;   // Force components

    Body(double m, double x_pos, double y_pos, double vel_x, double vel_y)
        : mass(m), x(x_pos), y(y_pos), vx(vel_x), vy(vel_y), fx(0.0), fy(0.0) {} //constructor

    void Force(const Body& other) { // force of another body on this body
        double dx = other.x - x;
        double dy = other.y - y;
        double dist = std::sqrt(sqr(dx) + sqr(dy));
        
        // To avoid division by zero
        dist = std::max(dist, NOT_ZERO);

        // Newton's law:
        double force = G * mass * other.mass / (dist * dist);
        
        // add these forces to the total exerted force
        fx += force * dx / dist;
        fy += force * dy / dist;
    }

    void updateVelocity(double dt) { // where dt is the timestep
        vx += dt * fx / mass;
        vy += dt * fy / mass;
    }

    void updatePosition(double dt) {
        x += dt * vx;
        y += dt * vy;
    }
    
};


class NBodySimulation {
private:
    std::vector<Body> bodies;
    double timeStep;
    double totalTime;
    double currentTime;
    int numSteps;

public:
    NBodySimulation(double dt, double total_time)
        : timeStep(dt), totalTime(total_time), currentTime(0.0) {
        numSteps = static_cast<int>(totalTime / timeStep);
    }
};