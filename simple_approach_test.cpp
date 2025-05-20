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
#include <iomanip>

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

    void newBody(const Body& body) {
        bodies.push_back(body);
    }

    // Simple approach
    void runSequential() {
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // Main loop
        for (int step = 1; step <= numSteps; ++step) {
            
            // reset the forces
            for (auto& body : bodies) {
                body.fx = 0.0;
                body.fy = 0.0;
            }
            
            // calculate forces between all pairs of bodies
            for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = 0; j < bodies.size(); ++j) {
                    if (i != j) {
                        bodies[i].Force(bodies[j]);
                    }
                }
            }
            
            // update velocities and positions
            for (auto& body : bodies) {
                body.updateVelocity(timeStep);
                body.updatePosition(timeStep);
            }
            
            // update current time
            currentTime += timeStep;
        
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Sequential simulation completed in " << duration << " ms" << std::endl;
    }
};

//////////////// to test
void createRandomSystem(NBodySimulation& sim, int numBodies) {
    // Central massive body (like a star)
    sim.newBody(Body(1.0e30, 0.0, 0.0, 0.0, 0.0));
    
    // Add random bodies
    for (int i = 1; i < numBodies; ++i) {
        double angle = (rand() % 1000) * 2.0 * M_PI / 1000.0;
        double distance = 1.0e11 + (rand() % 1000) * 1.0e9;
        double mass = 1.0e23 + (rand() % 1000) * 1.0e22;
        
        double x = distance * cos(angle);
        double y = distance * sin(angle);
        
        // Calculate circular orbit velocity
        double v = std::sqrt(G * 1.0e30 / distance);
        double vx = -v * sin(angle);
        double vy = v * cos(angle);
        
        sim.newBody(Body(mass, x, y, vx, vy));
    }
}
///////////////////////

int main() {
    // Set random seed
    srand(static_cast<unsigned int>(time(nullptr)));
    
    // Create simulation with time step and total time
    NBodySimulation simulation(1, 20); // 1 second time step, 20 total
    
    // Create a system of bodies 
    createRandomSystem(simulation, 100); 
    
    // Run sequential simulation
    simulation.runSequential();
    
    return 0;
}

/*
When compiling:  g++ -std=c++11 simple_approach_test.cpp -o test
When running: ./test

We obtain this output: "Sequential simulation completed in 7 ms"

*/