#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <iomanip>
#include <fstream>

#define G 6.67430e-11 // Gravity constant
#define NOT_ZERO 1e-9 // constant to avoid division by zero

double sqr(double x) {
    return x*x;
}

struct Vector {
    double data[3];
    Vector(double r = 0, double g = 0, double b = 0) {
        data[0] = r; data[1] = g; data[2] = b;
    }
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }
};

// Structure to represent a body in 2D space: we use a struct instead of a class as we want to make it publicly accessible in the whole code
class Body {
public:
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
        //double not_zero = NOT_ZERO;
        dist = std::max(dist, NOT_ZERO);

        // Newton's law:
        double force = G * mass * other.mass / (dist*dist);
        
        // compute the x and y components of the force
        double force_x = force * dx / dist;
        double force_y = force * dy / dist;

        //// By Newton's 3rd law Fij = -Fji
        // update current object's force
        this->fx += force_x;
        this->fy += force_y;

        // update the force of the object acting on it
        other.fx -= force_x;
        other.fy -= force_y;
    }

    void updateVelocity(double dt) { // where dt is the timestep
        this->vx += dt * fx / mass;
        this->vy += dt * fy / mass;
    }

    void updatePosition(double dt) {
        this->x += dt * vx;
        this->y += dt * vy;
    }
    
};

void save_frame_nbody(const std::vector<Body>& bodies, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 0); // Black background
    
    // Find bounds of all bodies to scale the view appropriately
    double minx = 1E30, miny = 1E30, maxx = -1E30, maxy = -1E30;
    for (const auto& body : bodies) {
        minx = std::min(minx, body.x);
        miny = std::min(miny, body.y);
        maxx = std::max(maxx, body.x);
        maxy = std::max(maxy, body.y);
    }
    
    // Add some padding (10% on each side)
    double rangex = maxx - minx;
    double rangey = maxy - miny;
    double padding = 0.1;
    minx -= rangex * padding;
    maxx += rangex * padding;
    miny -= rangey * padding;
    maxy += rangey * padding;
    
    // Update ranges after padding
    rangex = maxx - minx;
    rangey = maxy - miny;
    
    // Use the larger range to maintain aspect ratio
    double range = std::max(rangex, rangey);
    double centerx = (minx + maxx) / 2.0;
    double centery = (miny + maxy) / 2.0;
    
    // Draw each body
    for (size_t i = 0; i < bodies.size(); i++) {
        const Body& body = bodies[i];
        
        // Convert world coordinates to screen coordinates
        double screen_x = W * (body.x - centerx + range/2) / range;
        double screen_y = H * (body.y - centery + range/2) / range;
        
        // Determine body color based on mass (or you can use other criteria)
        Vector color;
        if (i == 0) {
            // Central massive body in yellow/orange
            color = Vector(255, 200, 0);
        } else {
            // Other bodies in different colors based on index
            double hue = (double)i / bodies.size() * 360.0;
            // Simple HSV to RGB conversion for variety
            int h_i = (int)(hue / 60) % 6;
            double f = hue / 60.0 - h_i;
            double p = 0, q = 1 - f, t = f;
            
            switch(h_i) {
                case 0: color = Vector(255, 255*t, 0); break;
                case 1: color = Vector(255*q, 255, 0); break;
                case 2: color = Vector(0, 255, 255*t); break;
                case 3: color = Vector(0, 255*q, 255); break;
                case 4: color = Vector(255*t, 0, 255); break;
                case 5: color = Vector(255, 0, 255*q); break;
            }
        }
        
        // Determine radius based on mass (logarithmic scale for better visualization)
        double base_radius = 3.0;
        double radius = base_radius;
        if (body.mass > 0) {
            radius = base_radius + 5.0 * log10(body.mass / 1e20); // Adjust scale as needed
            radius = std::max(2.0, std::min(radius, 20.0)); // Clamp between 2 and 20 pixels
        }
        
        // Draw filled circle
        int cx = (int)screen_x;
        int cy = (int)screen_y;
        int r = (int)radius;
        
        for (int dy = -r; dy <= r; dy++) {
            for (int dx = -r; dx <= r; dx++) {
                if (dx*dx + dy*dy <= r*r) {
                    int px = cx + dx;
                    int py = cy + dy;
                    
                    // Check bounds
                    if (px >= 0 && px < W && py >= 0 && py < H) {
                        int idx = ((H - py - 1) * W + px) * 3;
                        image[idx] = (unsigned char)color[0];
                        image[idx + 1] = (unsigned char)color[1];
                        image[idx + 2] = (unsigned char)color[2];
                    }
                }
            }
        }
        
        // Draw a smaller bright center for better visibility
        int inner_r = std::max(1, r/3);
        for (int dy = -inner_r; dy <= inner_r; dy++) {
            for (int dx = -inner_r; dx <= inner_r; dx++) {
                if (dx*dx + dy*dy <= inner_r*inner_r) {
                    int px = cx + dx;
                    int py = cy + dy;
                    
                    if (px >= 0 && px < W && py >= 0 && py < H) {
                        int idx = ((H - py - 1) * W + px) * 3;
                        image[idx] = 255;
                        image[idx + 1] = 255;
                        image[idx + 2] = 255;
                    }
                }
            }
        }
    }
    
    // Save the image
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}


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
        std::ofstream file("positions_sequential.csv");

        auto startTime = std::chrono::high_resolution_clock::now();
        save_frame_nbody(bodies, "images/sequential_frame_", 0);
        
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
                        std::cout << "Force on body " << i << ": (" << bodies[i].fx << ", " << bodies[i].fy << ")" << std::endl;
                    }
                }
            }
            
            // update velocities and positions
            for (auto& body : bodies) {
                body.updateVelocity(timeStep);
                body.updatePosition(timeStep);
                std::cout << "Body position: (" << body.x << ", " << body.y << ")" << std::endl;
                file << body.x << ',' << body.y << ',';
            }
            file << "\n";

            if (step % 10 == 0) {  // Save every 10th frame
                save_frame_nbody(bodies, "images/sequential_frame_", step / 10);
            }
            // update current time
            currentTime += timeStep;
        
        }
        save_frame_nbody(bodies, "images/sequential_frame_", numSteps / 10 + 1);
        file.close();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Sequential simulation completed in " << duration << " ms" << std::endl;
    }

    void updatePositionThread(std::vector<Body>::iterator begin, std::vector<Body>::iterator end, double dt) { 
        //std::ofstream file("position_parallel.csv");
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

        save_frame_nbody(bodies, "images/parallel_frame_", 0);
        
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

            if (step % 10 == 0) {  // Save every 10th frame
                save_frame_nbody(bodies, "images/parallel_frame_", step / 10);
            }

            currentTime += dt;
        }
        save_frame_nbody(bodies, "images/parallel_frame_", numSteps / 10 + 1);
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


        save_frame_nbody(bodies, "images/mutex_frame_", 0);


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

            if ((step + 1) % 10 == 0) {  // +1 because step starts from 0
                save_frame_nbody(bodies, "images/mutex_frame_", (step + 1) / 10);
            }

            currentTime += dt;
        }
        save_frame_nbody(bodies, "images/mutex_frame_", numSteps / 10 + 1);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "Parallel simulation with mutex completed in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << " ms" << std::endl;
    }


    void runParallelNoMutex(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();
        

        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return;
        if (Nthreads > length) Nthreads = length;

        save_frame_nbody(bodies, "images/nomutex_frame_", 0);

        for (int step = 0; step < numSteps; ++step) {
            std::vector<std::vector<std::pair<double, double>>> force_acc(length, std::vector<std::pair<double, double>>(Nthreads, {0.0, 0.0}));

            std::vector<std::thread> threads;
            auto computeForces = [&](size_t tid, size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    for (size_t j = 0; j < length; ++j) {
                        if (i == j) continue;
                        double dx = bodies[j].x - bodies[i].x;
                        double dy = bodies[j].y - bodies[i].y;
                        double dist = std::sqrt(sqr(dx) + sqr(dy));
                        dist = std::max(dist, NOT_ZERO);
                        double force = G * bodies[i].mass * bodies[j].mass / sqr(dist);
                        double fx = force * dx / dist;
                        double fy = force * dy / dist;
                        force_acc[i][tid].first += fx;
                        force_acc[i][tid].second += fy;
                    }
                }
            };

            size_t block_size = length / Nthreads;
            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end = (t == Nthreads - 1) ? length : start + block_size;
                threads.emplace_back(computeForces, t, start, end);
            }

            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            // Reduce forces
            for (size_t i = 0; i < length; ++i) {
                bodies[i].fx = bodies[i].fy = 0.0;
                for (size_t t = 0; t < Nthreads; ++t) {
                    bodies[i].fx += force_acc[i][t].first;
                    bodies[i].fy += force_acc[i][t].second;
                }
            }

            // Update positions and velocities
            threads.clear();
            auto updateFunc = [&](size_t start, size_t end) {
                updatePositionThread(bodies.begin() + start, bodies.begin() + end, dt);
            };

            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end = (t == Nthreads - 1) ? length : start + block_size;
                threads.emplace_back(updateFunc, start, end);
            }

            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            if ((step + 1) % 10 == 0) {  // +1 because step starts from 0
                save_frame_nbody(bodies, "images/nomutex_frame_", (step + 1) / 10);
            }

            currentTime += dt;
        }

        save_frame_nbody(bodies, "images/nomutex_frame_", numSteps / 10 + 1);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "Parallel simulation without mutex completed in "
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
        double v = std::sqrt(G * 1.0e30 / distance)*10000000;
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
    NBodySimulation simulation_sequential(1, 100); // 1 second time step, 20 total

    // Create a system of bodies 
    createRandomSystem(simulation_sequential, 10); 
    std::vector<Body> initial_bodies = simulation_sequential.getBodies();
    
    // Run sequential simulation
    simulation_sequential.runSequential();
    std::vector<Body> sequential_result = simulation_sequential.getBodies();

    NBodySimulation simulation_parallel(1, 100);
    simulation_parallel.setBodies(initial_bodies);
    simulation_parallel.runParallel(1.0, 2);
    std::vector<Body> parallel_result = simulation_parallel.getBodies();


    NBodySimulation simulation_parallel_mutex(1, 100);
    simulation_parallel_mutex.setBodies(initial_bodies);
    simulation_parallel_mutex.runParallelWithMutex(1.0, 2);
    std::vector<Body> mutex_result = simulation_parallel_mutex.getBodies();


    NBodySimulation simulation_parallel_nomutex(1, 100);
    simulation_parallel_nomutex.setBodies(initial_bodies);
    simulation_parallel_nomutex.runParallelNoMutex(1.0, 2);
    std::vector<Body> nomutex_result = simulation_parallel_nomutex.getBodies();

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
    if (compareResults(sequential_result, nomutex_result)) {
        std::cout << "OK." << std::endl;
    } else {
        std::cout << "No mutex version and sequential simulations grant different results" << std::endl;
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


// Alternative version that also draws trails
void save_frame_nbody_with_trails(const std::vector<Body>& bodies, 
                                 const std::vector<std::vector<std::pair<double, double>>>& trails,
                                 std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 0); // Black background
    
    // Find bounds (considering both current positions and trail history)
    double minx = 1E30, miny = 1E30, maxx = -1E30, maxy = -1E30;
    
    // Check current positions
    for (const auto& body : bodies) {
        minx = std::min(minx, body.x);
        miny = std::min(miny, body.y);
        maxx = std::max(maxx, body.x);
        maxy = std::max(maxy, body.y);
    }
    
    // Check trail positions
    for (const auto& trail : trails) {
        for (const auto& pos : trail) {
            minx = std::min(minx, pos.first);
            miny = std::min(miny, pos.second);
            maxx = std::max(maxx, pos.first);
            maxy = std::max(maxy, pos.second);
        }
    }
    
    // Add padding and calculate scaling
    double rangex = maxx - minx;
    double rangey = maxy - miny;
    double padding = 0.1;
    minx -= rangex * padding;
    maxx += rangex * padding;
    miny -= rangey * padding;
    maxy += rangey * padding;
    
    double range = std::max(maxx - minx, maxy - miny);
    double centerx = (minx + maxx) / 2.0;
    double centery = (miny + maxy) / 2.0;
    
    // Draw trails first (so they appear behind the bodies)
    for (size_t i = 0; i < trails.size() && i < bodies.size(); i++) {
        const auto& trail = trails[i];
        
        // Color for trail (dimmer version of body color)
        Vector trail_color;
        if (i == 0) {
            trail_color = Vector(100, 80, 0); // Dim yellow
        } else {
            double hue = (double)i / bodies.size() * 360.0;
            int h_i = (int)(hue / 60) % 6;
            double f = hue / 60.0 - h_i;
            double p = 0, q = 1 - f, t = f;
            
            switch(h_i) {
                case 0: trail_color = Vector(100, 100*t, 0); break;
                case 1: trail_color = Vector(100*q, 100, 0); break;
                case 2: trail_color = Vector(0, 100, 100*t); break;
                case 3: trail_color = Vector(0, 100*q, 100); break;
                case 4: trail_color = Vector(100*t, 0, 100); break;
                case 5: trail_color = Vector(100, 0, 100*q); break;
            }
        }
        
        // Draw trail points
        for (const auto& pos : trail) {
            double screen_x = W * (pos.first - centerx + range/2) / range;
            double screen_y = H * (pos.second - centery + range/2) / range;
            
            int px = (int)screen_x;
            int py = (int)screen_y;
            
            if (px >= 0 && px < W && py >= 0 && py < H) {
                int idx = ((H - py - 1) * W + px) * 3;
                image[idx] = (unsigned char)trail_color[0];
                image[idx + 1] = (unsigned char)trail_color[1];
                image[idx + 2] = (unsigned char)trail_color[2];
            }
        }
    }
    
    // Draw bodies (same as previous function)
    for (size_t i = 0; i < bodies.size(); i++) {
        const Body& body = bodies[i];
        
        double screen_x = W * (body.x - centerx + range/2) / range;
        double screen_y = H * (body.y - centery + range/2) / range;
        
        Vector color;
        if (i == 0) {
            color = Vector(255, 200, 0);
        } else {
            double hue = (double)i / bodies.size() * 360.0;
            int h_i = (int)(hue / 60) % 6;
            double f = hue / 60.0 - h_i;
            double p = 0, q = 1 - f, t = f;
            
            switch(h_i) {
                case 0: color = Vector(255, 255*t, 0); break;
                case 1: color = Vector(255*q, 255, 0); break;
                case 2: color = Vector(0, 255, 255*t); break;
                case 3: color = Vector(0, 255*q, 255); break;
                case 4: color = Vector(255*t, 0, 255); break;
                case 5: color = Vector(255, 0, 255*q); break;
            }
        }
        
        double base_radius = 3.0;
        double radius = base_radius;
        if (body.mass > 0) {
            radius = base_radius + 5.0 * log10(body.mass / 1e20);
            radius = std::max(2.0, std::min(radius, 20.0));
        }
        
        int cx = (int)screen_x;
        int cy = (int)screen_y;
        int r = (int)radius;
        
        // Draw filled circle
        for (int dy = -r; dy <= r; dy++) {
            for (int dx = -r; dx <= r; dx++) {
                if (dx*dx + dy*dy <= r*r) {
                    int px = cx + dx;
                    int py = cy + dy;
                    
                    if (px >= 0 && px < W && py >= 0 && py < H) {
                        int idx = ((H - py - 1) * W + px) * 3;
                        image[idx] = (unsigned char)color[0];
                        image[idx + 1] = (unsigned char)color[1];
                        image[idx + 2] = (unsigned char)color[2];
                    }
                }
            }
        }
        
        // Bright center
        int inner_r = std::max(1, r/3);
        for (int dy = -inner_r; dy <= inner_r; dy++) {
            for (int dx = -inner_r; dx <= inner_r; dx++) {
                if (dx*dx + dy*dy <= inner_r*inner_r) {
                    int px = cx + dx;
                    int py = cy + dy;
                    
                    if (px >= 0 && px < W && py >= 0 && py < H) {
                        int idx = ((H - py - 1) * W + px) * 3;
                        image[idx] = 255;
                        image[idx + 1] = 255;
                        image[idx + 2] = 255;
                    }
                }
            }
        }
    }
    
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}
