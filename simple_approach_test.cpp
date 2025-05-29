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
        double force = G * mass * other.mass / sqr(dist);
        
        // add these forces to the total exerted force
        fx += force * dx / dist;
        fy += force * dy / dist;
    }

    void OptimizedForce(Body& other) {
        // Directly computes Fij and Fji
        double dx = other.x - x;
        double dy = other.y - y;
        double dist = std::sqrt(sqr(dx) + sqr(dy));
        
        // To avoid division by zero
        dist = std::max(dist, NOT_ZERO);

        // Newton's law:
        double force = G * mass * other.mass / sqr(dist);
        
        // compute the x and y components of the force
        double force_x = force * dx / dist;
        double force_y = force * dy / dist;

        //// By Newton's 3rd law Fij = -Fji
        // update current object's force
        fx += force_x;
        fy += force_y;

        // update the force of the object acting on it
        other.fx -= force_x;
        other.fy -= force_y;
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

    const std::vector<Body>& getBodies() { 
        return bodies; 
    }

    void setBodies(const std::vector<Body>& newBodies) { 
        bodies = newBodies; 
    }

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

    void updatePositionThread(std::vector<Body>::iterator begin, std::vector<Body>::iterator end, double dt) { 
        std::vector<Body>::iterator it = begin;
        while (it != end) {
            it->updateVelocity(dt);
            it->updatePosition(dt);
            ++it;
        }
}
    void runParallel(double dt, size_t Nthreads) {
        // Start timer
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();

        if (length == 0){
            return;
        } 
        if (Nthreads == 0) {
            Nthreads = 1;
        }
        
        for (int step = 1; step <= numSteps; ++step) {
            // Compute forces sequentially
            for (size_t i = 0; i < bodies.size(); ++i) {
                bodies[i].fx = 0.0;
                bodies[i].fy = 0.0;
            }
        
            for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = i+1; j < bodies.size(); ++j) { // originally we started at size_t j = 0
                    if (i != j) {
                        bodies[i].OptimizedForce(bodies[j]); // originally we used Force(bodies[j])
                    }
                }
            }

            // Update position and velocity in parallel
            std::vector<std::thread> threads(Nthreads - 1);
            std::vector<Body>::iterator block_start = bodies.begin();

            size_t block_size = length / Nthreads;

            for (size_t i = 0; i < Nthreads - 1; ++i) {
                std::vector<Body>::iterator block_end = block_start + block_size;
                threads[i] = std::thread(&NBodySimulation::updatePositionThread, this, block_start, block_end, dt); // chatGPT helped me figure out using "this"
                block_start = block_end;
            }

            updatePositionThread(block_start, bodies.end(), dt);

            //  Join Threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            currentTime += dt;
        }

        // End timer
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " ms" << std::endl;
    }

    void runParallelWithMutex(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return; // no bodies or threads
        if (Nthreads > length) Nthreads = length;

        std::vector<std::mutex> force_mutexes(length);

        for (int step = 0; step < numSteps; ++step) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                bodies[i].fx = 0.0;
                bodies[i].fy = 0.0;
            }

            std::vector<std::thread> threads;
            auto computeForces = [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    for (size_t j = 0; j < bodies.size(); ++j) {
                        if (i == j) continue;
                        double dx = bodies[j].x - bodies[i].x;
                        double dy = bodies[j].y - bodies[i].y;
                        double dist = std::sqrt(sqr(dx) + sqr(dy));
                        dist = std::max(dist, NOT_ZERO);
                        double force = G * bodies[i].mass * bodies[j].mass / sqr(dist);
                        double fx = force * dx / dist;
                        double fy = force * dy / dist;

                        std::lock_guard<std::mutex> lock(force_mutexes[i]);
                        bodies[i].fx += fx;
                        bodies[i].fy += fy;
                    }
                }
            };

            size_t block_size = length / Nthreads;
            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end = (t == Nthreads - 1) ? length : start + block_size;
                threads.emplace_back(computeForces, start, end);
            }

            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            // update positions and velocities
            std::vector<std::thread> updateThreads;
            auto updateFunc = [&](size_t start, size_t end) {
                updatePositionThread(bodies.begin() + start, bodies.begin() + end, dt);
            };

            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end = (t == Nthreads - 1) ? length : start + block_size;
                updateThreads.emplace_back(updateFunc, start, end);
            }

            for (size_t i = 0; i < updateThreads.size(); ++i) {
                updateThreads[i].join();
            }

            currentTime += dt;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "Parallel simulation with mutex completed in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << " ms" << std::endl;
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
    NBodySimulation simulation_sequential(1, 20); // 1 second time step, 20 total

    // Create a system of bodies 
    createRandomSystem(simulation_sequential, 1000); 
    std::vector<Body> initial_bodies = simulation_sequential.getBodies();
    
    // Run sequential simulation
    simulation_sequential.runSequential();
    std::vector<Body> sequential_result = simulation_sequential.getBodies();

    NBodySimulation simulation_parallel(1, 20);
    simulation_parallel.setBodies(initial_bodies);
    simulation_parallel.runParallel(1.0, 4);
    std::vector<Body> parallel_result = simulation_parallel.getBodies();


    NBodySimulation simulation_parallel_mutex(1, 20);
    simulation_parallel_mutex.setBodies(initial_bodies);
    simulation_parallel_mutex.runParallelWithMutex(1.0, 4);
    std::vector<Body> mutex_result = simulation_parallel_mutex.getBodies();

    //Test if both simulations grant the same result (ChatGPT helped me debug the code by adding the comparison of sizes before checking values and added the 'break')
    bool same_results = true;

    if (sequential_result.size() != parallel_result.size()) {
        same_results = false;
    } 
    
    else {
        for (size_t i=0; i < sequential_result.size(); ++i) {
            double dx = std::abs(sequential_result[i].x - parallel_result[i].x);
            double dy = std::abs(sequential_result[i].y - parallel_result[i].y);
            double dvx = std::abs(sequential_result[i].vx - parallel_result[i].vx);
            double dvy = std::abs(sequential_result[i].vy - parallel_result[i].vy);

            if (dx > 1e-6 || dy > 1e-6 || dvx > 1e-6 || dvy > 1e-6) {
                same_results = false;
                break;
            }
        }
    }

    if (same_results) {
        std::cout << "OK." << std::endl;
    } else {
        std::cout << "Parallel and sequential simulations grant different results" << std::endl;
    }



    ////////////// 
    // Compare results of parallel and mutex 
    auto compareResults = [](const std::vector<Body>& a, const std::vector<Body>& b) {
        if (a.size() != b.size()) return false;

        for (size_t i = 0; i < a.size(); ++i) {
            double dx = std::abs(a[i].x - b[i].x);
            double dy = std::abs(a[i].y - b[i].y);
            double dvx = std::abs(a[i].vx - b[i].vx);
            double dvy = std::abs(a[i].vy - b[i].vy);

            if (dx > 1e-6 || dy > 1e-6 || dvx > 1e-6 || dvy > 1e-6) {
                return false;
            }
        }
        return true;
    };
    if (compareResults(sequential_result, mutex_result)) {
        std::cout << "OK." << std::endl;
    } else {
        std::cout << "Mutex and sequential simulations grant different results" << std::endl;
    }
    
    return 0;
}

/*
When compiling:  g++ -std=c++11 simple_approach_test.cpp -o test
When running: ./test

1) When we only had runSequential
Sequential simulation completed in 7 ms

2) When we added parallelization
Sequential simulation completed in 7 ms
Parallel simulation completed in 13 ms

3) When we avoided computing Fij twice
Sequential simulation completed in 9 ms
Parallel simulation completed in 7 ms

*/