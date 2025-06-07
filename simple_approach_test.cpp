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
#include <fstream>

////////////////////////////////////////////////////////////////////////////////////////
//////////////// DEFINING CONSTANTS AND USEFUL FUNCTIONS  //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// Constants for all code
#define G 6.67430e-11 // Gravity constant
#define NOT_ZERO 1e-9 // constant to avoid division by zero

// Constants for Barnes-Hutt Algorithm
const int MAX_BODIES = 10000; 
const int MAX_NODES = 10 * MAX_BODIES;
const double THETA = 0.9;                 // Barnes-Hut threshold
const double SOFTENING = 1e-3;            // Suggested as an improvement by chatgpt, avoid infinity when two bodies get very close to each other
const double MIN_SUBDIVIDE = 1e-1;        // This can be adapted to the system we study (should be smaller than the expected distances between bodies, 100 seems adapted for the solar system)

// Square function
double sqr(double x) {
    return x*x;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////// DEFINING BODY CLASS TO MANIPULATE THE OBJECT BODIES ///////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// Body class to store data for each of the N bodies
class Body {
public:
    double mass;     // Mass of the body
    double x, y;     // Position coordinates
    double vx, vy;   // Velocity components
    double fx, fy;   // Force components

    Body(double m, double x_pos, double y_pos, double vel_x, double vel_y)
        : mass(m), x(x_pos), y(y_pos), vx(vel_x), vy(vel_y), fx(0.0), fy(0.0) {} //constructor

    void Force(const Body& other) { // force of another body on this body
        double dx = other.x - x;
        double dy = other.y - y;
        double dist = std::sqrt(sqr(dx) + sqr(dy));
        
        dist = std::max(dist, NOT_ZERO); // avoid division by zero

        // Newton's law:
        double force = G * mass * other.mass / sqr(dist);
        
        // add these forces to the total exerted force
        fx += force * dx / dist;
        fy += force * dy / dist;
    }

    void OptimizedForce(Body& other) {
        // directly computes Fij and Fji
        double dx = other.x - x;
        double dy = other.y - y;
        double dist = std::sqrt(sqr(dx) + sqr(dy));
        
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

////////////////////////////////////////////////////////////////////////////////////////
//////////// DEFINING QUADTREE STRUCTURE FOR THE BARNES-HUTT ALGORITHM ////////////////
////////////////////////////////////////////////////////////////////////////////////////

// Nodes of the quadtree for the Barnes-Hut algorithm
struct QuadNode {
    bool hasBody = false;
    int children[4] = {-1, -1, -1, -1};
    int bodyIndex = -1;

    double centerX, centerY; // center or the region the node is in
    double halfSize; // half of a side of the square it is in
    double mass = 0.0; // total mass of the region
    double comX = 0.0; // center of mass of the region
    double comY = 0.0;
};

// Quadtree for the Barnes-Hut algorithm
class Quadtree {
public:
    Quadtree() : nodeCount(0) {} // By default the Quatree has 0 nodes

    QuadNode quadtree[MAX_NODES]; // Quadtrees are implemented as arrays with a set length
    int nodeCount; // Counter to ensure each node is in a different cell of the quadtree array

    // Reset to empty tree
    void reset() { 
        nodeCount = 0; 
    }

    // Create a new node
    int createNode(double cx, double cy, double hs) {
        // Check if there is still room for an extra node in the tree array
        if (nodeCount >= MAX_NODES) {
            std::cerr << "Error: Maximum number of nodes exceeded\n";
            std::exit(1); // interrupt the program directly
        }
        // Initialize information of the region of space to which the node belongs
        int id = nodeCount++;
        quadtree[id] = QuadNode();
        quadtree[id].centerX = cx;
        quadtree[id].centerY = cy;
        quadtree[id].halfSize = hs;
        return id;
    }

    // Determine which quadrant a body b belongs to (value between 0 and 3)
    int getQuadrant(const QuadNode& node, const Body& b) {
        int quad = 0;
        if (b.x > node.centerX) quad += 1;
        if (b.y > node.centerY) quad += 2;
        return quad;
    }

    // Insert a body 'b' of index 'body_index' into the tree rooted at node 'x' of index 'node_index'
    void insert(int node_index, int body_index, std::vector<Body>& bodies) {
        QuadNode& node_x = quadtree[node_index];
        Body& body_b = bodies[body_index];

        // 0. Check that you're not subdiving too much - added after encountering infinite loops due to poorly chosen domainSize
        if (node_x.halfSize <= MIN_SUBDIVIDE) {
            std::cerr << "Error: The space has been subdivided below 1e3, check you have chosen the write domainSize. Otherwise adjust MIN_SUBDIVIDE\n";
            return;
        }
        // 1. If node x does not contain a body, 
        if (!node_x.hasBody && node_x.bodyIndex == -1 && node_x.children[0] == -1) { // double check there is no body + check that it is not an internal node
            // put the new body b here.
            node_x.bodyIndex = body_index;
            node_x.hasBody = true;
            node_x.mass = body_b.mass;
            node_x.comX = body_b.x;
            node_x.comY = body_b.y;
            return;
        } 
        else {
            // 3. If node x is an external node, say containing a body named c.   
            // Since b and c may still end up in the same quadrant, there may be several subdivisions during a single insertion. 
            // Finally, update the center-of-mass and total mass of x.
            if (node_x.hasBody) { // then it is necesseraly an external node containing a body
                // then there are two bodies b and c in the same region.
                int body_c_index = node_x.bodyIndex;
                Body& body_c = bodies[body_c_index];
                // we take body_c out of node_x
                node_x.bodyIndex = -1;
                node_x.hasBody = false;

                // Subdivide the region further by creating four children. chatGPT coded this part.
                for (int i = 0; i < 4; ++i) {
                    double offsetX = (i % 2 == 0 ? -0.5 : 0.5) * node_x.halfSize;
                    double offsetY = (i < 2 ? -0.5 : 0.5) * node_x.halfSize;
                    node_x.children[i] = createNode(
                        node_x.centerX + offsetX,
                        node_x.centerY + offsetY,
                        node_x.halfSize / 2
                    );
                }

                // Then, recursively insert both b and c into the appropriate quadrant(s).
                int oldQuad = getQuadrant(node_x, body_c); // find the quadrant in which c was in
                insert(node_x.children[oldQuad], body_c_index, bodies); // insert body c into a child quadrant
            }
            // 3. and 2. If node x is an internal node, recursively insert the body b in the appropriate quadrant.
            int quad = getQuadrant(node_x, body_b);
            insert(node_x.children[quad], body_index, bodies);
        }
        // Update total mass and center of mass
        double totalMass = 0.0;
        double weightedX = 0.0;
        double weightedY = 0.0;
        for (int i=0; i<4; ++i) {
            int childID = node_x.children[i];
            if (childID == -1) continue;
            QuadNode& child = quadtree[childID];
            totalMass += child.mass;
            weightedX += child.comX*child.mass;
            weightedY += child.comY*child.mass;
        }
        node_x.mass = totalMass;
        node_x.comX = weightedX/totalMass;
        node_x.comY = weightedY/totalMass;
    }

    // Calculate the net force acting on body b sequentially
    void computeForce(int nodeID, Body& b, std::vector<Body>& bodies) {
        // Uninitialized node
        if (nodeID == -1) {
            return;
        }

        QuadNode& node = quadtree[nodeID];
        // Node with 0 mass
        if (node.mass == 0) return;

        double dx = node.comX - b.x;
        double dy = node.comY - b.y;
        double dist = sqrt(dx * dx + dy * dy + SOFTENING * SOFTENING);

        // If the current node is an external node and is the body b
        if (node.hasBody && node.bodyIndex == (&b - &bodies[0])) {
            return;
        }

        // If the current node is an external node with a body different from b and s/d < θ
        if (node.hasBody || node.halfSize / dist < THETA) {
            // calculate the force it exerts on body b
            double F = G * b.mass * node.mass / (dist * dist + SOFTENING * SOFTENING);
            // add this amount to b’s net force
            b.fx += F * dx / dist;
            b.fy += F * dy / dist;

        // If the current node is an external node with a body different from b and s/d >= θ
        } else {
            // run the procedure recursively on each of the current node’s children.
            for (int i = 0; i < 4; ++i) {
                computeForce(node.children[i], b, bodies);
            }
        }
    }

    // Calculate the net force acting on body b thread version
    void computeForceThread(std::vector<Body>::iterator begin, std::vector<Body>::iterator end, int root, std::vector<Body>& bodies) {
        for (auto it=begin; it!=end; ++it) {
            it->fx = 0.0;
            it->fy = 0.0;
            computeForce(root, *it, bodies);
        }
    }

    // Update position and velocity of the bodies sequencially
    void update_pos_vel(std::vector<Body>& bodies, double dt) {
        for (Body& b : bodies) {
            b.updatePosition(dt);
            b.updateVelocity(dt);
        }
    }

    // Update position and velocity of the bodies threaded version
    void update_pos_vel_Thread(std::vector<Body>::iterator begin, std::vector<Body>::iterator end, std::vector<Body>& bodies, double dt) {
        for (auto it=begin; it!=end; ++it) {
            it->updatePosition(dt);
            it->updateVelocity(dt);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// CLASS TO RUN ALL SIMULATIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class NBodySimulation {
private:
    std::vector<Body> bodies; // store all bodies
    double timeStep; 
    double totalTime; // total simulation time
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

    void newBody(const Body& body) { // add to bodies
        bodies.push_back(body);
    }

    ////// BEGINNING OF DIFFERENT SIMULATION METHODS ///////////////////////////////////////////////////////////////////
    
    // Simple sequential approach without any form of parallelization ==================================================
    long runSequential() {
        std::ofstream file("positions_sequential.csv"); // store the bodies' positions at each timestep for plotting
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            return 0;
        }

        for (size_t i = 0; i < bodies.size(); ++i) {
            file << bodies[i].mass;
            if (i != bodies.size() - 1) file << ',';
        }
        file << '\n';

        auto startTime = std::chrono::high_resolution_clock::now();
        
        // main loop
        for (int step = 1; step <= numSteps; ++step) {
            
            // reset the forces
            for (auto& body : bodies) {
                body.fx = 0.0;
                body.fy = 0.0;
            }
            
            // calculate forces between all pairs of bodies
            /*for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = 0; j < bodies.size(); ++j) {
                    if (i != j) {
                        bodies[i].Force(bodies[j]);
                        std::cout << "Force on body " << i << ": (" << bodies[i].fx << ", " << bodies[i].fy << ")" << std::endl;
                    }
                }
            }*/
            for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = i + 1; j < bodies.size(); ++j) {
                    bodies[i].OptimizedForce(bodies[j]);
                }
            }

            
            // update velocities and positions
            for (auto& body : bodies) {
                body.updateVelocity(timeStep);
                body.updatePosition(timeStep);
                //std::cout << "Body position: (" << body.x << ", " << body.y << ")" << std::endl;
                file << body.x << ',' << body.y << ',';
            }
            file << "\n";

            // update current time
            currentTime += timeStep;
        
        }
        file.close();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Sequential simulation completed in " << duration << " ms" << std::endl;
        return duration;
    }

    // First parallelization: updating the position and velocity in parallel =================================
    void updatePositionThread(std::vector<Body>::iterator begin, std::vector<Body>::iterator end, double dt) { 
        //std::ofstream file("position_parallel.csv");
        std::vector<Body>::iterator it = begin;
        while (it != end) {
            it->updateVelocity(dt);
            it->updatePosition(dt);
            ++it;
        }
    }
    long runParallel(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();

        if (length == 0){
            return 0;
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
            // update position and velocity in parallel
            std::vector<std::thread> threads(Nthreads - 1);
            std::vector<Body>::iterator block_start = bodies.begin();
            size_t block_size = length / Nthreads;
            // Update position and velocity in parallel
            for (size_t i = 0; i < Nthreads - 1; ++i) {
                std::vector<Body>::iterator block_end = block_start + block_size;
                threads[i] = std::thread(&NBodySimulation::updatePositionThread, this, block_start, block_end, dt); // chatGPT helped me figure out using "this"
                block_start = block_end;
            }
            // Last thread to update posotion and velocity
            updatePositionThread(block_start, bodies.end(), dt);
            // Join threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }
            currentTime += dt;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " ms" << std::endl;
        return duration;
    }

    long runParallelWithMutex(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return 0; // no bodies or threads
        if (Nthreads > length) Nthreads = length;

        std::vector<std::mutex> force_mutexes(length); // mutex for each body to not have issues with the computation of its forces

        for (int step = 0; step < numSteps; ++step) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                bodies[i].fx = 0.0;
                bodies[i].fy = 0.0;
            }

            // compute forces in parallel
            std::vector<std::thread> threads;
            auto computeForces = [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    for (size_t j = 0; j < bodies.size(); ++j) {
                        if (i == j) continue; // do not compute force with itself
                        double dx = bodies[j].x - bodies[i].x;
                        double dy = bodies[j].y - bodies[i].y;
                        double dist = std::sqrt(sqr(dx) + sqr(dy));
                        dist = std::max(dist, NOT_ZERO);
                        double force = G * bodies[i].mass * bodies[j].mass / sqr(dist);
                        double fx = force * dx / dist;
                        double fy = force * dy / dist;

                        std::lock_guard<std::mutex> lock(force_mutexes[i]); // update with mutex to protect forces
                        bodies[i].fx += fx;
                        bodies[i].fy += fy;
                    }
                }
            };

            size_t block_size = length / Nthreads;
            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
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
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                updateThreads.emplace_back(updateFunc, start, end);
            }

            for (size_t i = 0; i < updateThreads.size(); ++i) {
                updateThreads[i].join();
            }

            currentTime += dt;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation with mutex completed in " << duration << " ms" << std::endl;
        return duration;
    }

    long runParallelNoMutex(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return 0;
        if (Nthreads > length) Nthreads = length;

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
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(computeForces, t, start, end);
            }

            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            // reduce forces
            for (size_t i = 0; i < length; ++i) {
                bodies[i].fx = bodies[i].fy = 0.0;
                for (size_t t = 0; t < Nthreads; ++t) {
                    bodies[i].fx += force_acc[i][t].first;
                    bodies[i].fy += force_acc[i][t].second;
                }
            }

            // update positions and velocities
            threads.clear();
            auto updateFunc = [&](size_t start, size_t end) {
                updatePositionThread(bodies.begin() + start, bodies.begin() + end, dt);
            };

            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(updateFunc, start, end);
            }

            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            currentTime += dt;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation without mutex completed in " << duration << " ms" << std::endl;
        return duration;
    }

    // Regular Barnes-Hut algorithm ============================================================================
    long runBarnesHutt(double domainSize, double dt) {
        Quadtree tree;

        auto startTime = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < numSteps; ++step) {
            tree.reset();
            int root = tree.createNode(0.0, 0.0, domainSize / 2.0);

            for (int i = 0; i < bodies.size(); ++i) {
                tree.insert(root, i, bodies);
            }

            for (Body& b : bodies) {
                b.fx = b.fy = 0.0;
                tree.computeForce(root, b, bodies);
            }
            tree.update_pos_vel(bodies, dt);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Barnes Hutt simulation completed in " << duration << " ms" << std::endl;
        return duration;
    }

    // Barnes-Hut algorithm: with parallelized force computation and position and velocity update  =================================
    void runBarnesHutParallel_Force_Position(double domainSize, double dt, size_t Nthreads) {
        Quadtree tree;
        size_t length = bodies.size();

        for (int step = 0; step < numSteps; ++step) {
            // Building tree sequentially
            tree.reset();
            int root = tree.createNode(0.0, 0.0, domainSize / 2.0);

            for (int i = 0; i < bodies.size(); ++i) {
                tree.insert(root, i, bodies);
            }
            
            // Initialize threads
            std::vector<std::thread> threads(Nthreads - 1);
            std::vector<Body>::iterator block_start = bodies.begin();

            if (length == 0){
                return;
            } 
            if (Nthreads == 0) {
                Nthreads = 1;
            }
            size_t block_size = length / Nthreads;

            // Updating forces in parallel
            for (size_t i = 0; i < Nthreads - 1; ++i) {
                std::vector<Body>::iterator block_end = block_start + block_size;
                threads[i] = std::thread(&Quadtree::computeForceThread, &tree, block_start, block_end, root, std::ref(bodies));
                block_start = block_end;
            }
            tree.computeForceThread(block_start, bodies.end(), root, bodies);
            
            // Join threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }
            
            // Clear the threads for reuse
            threads.clear(); // This was suggested to me by chatgpt as an improvement to reuse threads
            block_start = bodies.begin();

            // Update position and velocity in parallel
            for (size_t i = 0; i < Nthreads - 1; ++i) {
                std::vector<Body>::iterator block_end = block_start + block_size;
                threads[i] = std::thread(&Quadtree::update_pos_vel_Thread, &tree, block_start, block_end, std::ref(bodies), dt);
                block_start = block_end;
            }
            tree.update_pos_vel_Thread(block_start, bodies.end(), bodies, dt);
            
            // Join threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            // Update time
            currentTime += dt;
        }
    }
    
    // Barnes-Hut algorithm: as above + parallelized tree construction =================================
    void buildTreeThread(Quadtree& tree_quadranti, int& root, std::vector<Body>& bodies, double centerX, double centerY, double halfSize) {
        root = tree_quadranti.createNode(centerX, centerY, halfSize);
        for (size_t j = 0; j < bodies.size(); ++j) {
            tree_quadranti.insert(root, j, bodies);
        }
    }

    void buildTreeParallel(Quadtree& tree, int root, std::vector<Body>& bodies, double centerX, double centerY, double halfSize, size_t Nthreads) {
        Nthreads = std::min(Nthreads, size_t(4)); // we want at most 4 threads for 4 quadrants

        // Sort the N bodies by quadrants of the space
        std::vector<Body> bodies_in_quad[4]; // array of vector of bodies, each vector represents the bodies in a given quadrant of space
        for (Body& b : bodies) {
            int quad=0;
            if (b.x > centerX) quad+=1;
            if (b.y > centerY) quad+=2;
            bodies_in_quad[quad].push_back(b);
        }

        // Create arrays for the 4 quadtrees, their 4 roots and the threads
        std::array<Quadtree, 4> subTrees;
        std::array<int, 4> subtreeRoots= { -1, -1, -1, -1 };
        std::vector<std::thread> threads(Nthreads - 1);
        
        // Build 4 trees corresponding to the 4 quadrants of the space
        for (size_t i = 0; i < Nthreads - 1; ++i) {
            // Chatgpt coded the following
            double offsetX = (i % 2 == 0 ? -0.5 : 0.5) * halfSize;
            double offsetY = (i < 2 ? -0.5 : 0.5) * halfSize;
            double qCenterX = centerX + offsetX;
            double qCenterY = centerY + offsetY;
            threads[i] = std::thread(&NBodySimulation::buildTreeThread, this, std::ref(subTrees[i]), std::ref(subtreeRoots[i]), std::ref(bodies_in_quad[i]), qCenterX, qCenterY, halfSize / 2);
        }

        // Run the last thread to build the 4 trees
        if ((Nthreads-1)<4) {
            size_t thread = Nthreads-1;
            // Chatgpt coded the following
            double offsetX = (thread % 2 == 0 ? -0.5 : 0.5) * halfSize;
            double offsetY = (thread < 2 ? -0.5 : 0.5) * halfSize;
            double qCenterX = centerX + offsetX;
            double qCenterY = centerY + offsetY;
            this->buildTreeThread(subTrees[thread], subtreeRoots[thread], bodies_in_quad[thread],qCenterX, qCenterY, halfSize/2);
        }

        // Join threads
        for (size_t i = 0; i < threads.size(); ++i) {
            threads[i].join();
        }

        // Insert the 4 trees as children of the final tree, update the mass and com of the root
        QuadNode& root_node = tree.quadtree[root];
        root_node.mass = 0;
        root_node.comX = 0;
        root_node.comY = 0;
        for (size_t quadrant_i = 0; quadrant_i < 4; ++quadrant_i) {
            // Check if there are bodies in this quadrant
            if (subtreeRoots[quadrant_i] == -1) {
                root_node.children[quadrant_i] = -1;
                continue;
            }
            // Update indices of nodes to ensure unicity of the indices in the final tree
            int minNode = tree.nodeCount;
            for (int n = 0; n < subTrees[quadrant_i].nodeCount; ++n) {
                tree.quadtree[tree.nodeCount++] = subTrees[quadrant_i].quadtree[n];
            }
            // Link this subtree as a child of the root
            int child_i = subtreeRoots[quadrant_i] + minNode;
            root_node.children[quadrant_i] = child_i;
            // Update massand com of theroot
            QuadNode& child = tree.quadtree[child_i];
            root_node.mass += child.mass;
            root_node.comX += child.comX * child.mass;
            root_node.comY += child.comY * child.mass;
        }
        if (root_node.mass > 0) {
            root_node.comX /= root_node.mass;
            root_node.comY /= root_node.mass;
        }
    }

    void runBarnesHutParallel(double domainSize, double dt, size_t Nthreads) {
        Quadtree tree;
        size_t length = bodies.size();

        // Create threads
        if (Nthreads == 0) {
            Nthreads = 1;
        }
        if (length == 0){
            return;
        } 
        
        for (int step = 0; step < numSteps; ++step) {
            // Building tree in parallel
            tree.reset();
            int root = tree.createNode(0.0, 0.0, domainSize / 2.0);

            buildTreeParallel(tree, root, std::ref(bodies), 0.0, 0.0, domainSize / 2.0, Nthreads);
            
            // Prepare threads
            std::vector<std::thread> threads(Nthreads - 1);
            std::vector<Body>::iterator block_start = bodies.begin();
            
            size_t block_size = length / Nthreads;

            // Compute forces in parallel
            for (size_t i = 0; i < Nthreads - 1; ++i) {
                std::vector<Body>::iterator block_end = block_start + block_size;
                threads[i] = std::thread(&Quadtree::computeForceThread, &tree, block_start, block_end, root, std::ref(bodies));
                block_start = block_end;
            }
            tree.computeForceThread(block_start, bodies.end(), root, bodies);
            
            // Join threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }
            
            // Clear threads for reuse
            threads.clear(); // This was suggested to me by chatgpt as an improvement to reuse threads
            block_start = bodies.begin();

            // Update positions and velocities in parallel
            for (size_t i = 0; i < Nthreads - 1; ++i) {
                auto block_end = block_start + block_size;
                threads.emplace_back(&Quadtree::update_pos_vel_Thread, &tree, block_start, block_end, std::ref(bodies), dt);
                block_start = block_end;
            }
            tree.update_pos_vel_Thread(block_start, bodies.end(), bodies, dt);
            
            // Join threads
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }
            // Update time
            currentTime += dt;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////// FOR ACCURACY- we need a history of each timestep //////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Body>> runSequentialWithTrace() {
        std::vector<std::vector<Body>> history;

        auto startTime = std::chrono::high_resolution_clock::now();
        
        // main loop
        for (int step = 1; step <= numSteps; ++step) {
            
            // reset the forces
            for (auto& body : bodies) {
                body.fx = 0.0;
                body.fy = 0.0;
            }
            
            for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = i + 1; j < bodies.size(); ++j) {
                    bodies[i].OptimizedForce(bodies[j]);
                }
            }

            
            // update velocities and positions
            for (auto& body : bodies) {
                body.updateVelocity(timeStep);
                body.updatePosition(timeStep);
                //std::cout << "Body position: (" << body.x << ", " << body.y << ")" << std::endl;
            }

            history.push_back(bodies);

            // update current time
            currentTime += timeStep;
        
        }


        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Sequential simulation completed in " << duration << " ms" << std::endl;
        return history;
    }

    
    std::vector<std::vector<Body>> runParallelWithTrace(double dt, size_t Nthreads) {
        std::vector<std::vector<Body>> history;
        auto startTime = std::chrono::high_resolution_clock::now();
        
        size_t length = bodies.size();
        if (length == 0) return history;
        if (Nthreads == 0) Nthreads = 1;

        for (int step = 1; step <= numSteps; ++step) {
            for (auto& b : bodies) b.fx = b.fy = 0.0;

            for (size_t i = 0; i < bodies.size(); ++i) {
                for (size_t j = i+1; j < bodies.size(); ++j) {
                    bodies[i].OptimizedForce(bodies[j]);
                }
            }

            std::vector<std::thread> threads(Nthreads - 1);
            auto block_start = bodies.begin();
            size_t block_size = length / Nthreads;

            for (size_t i = 0; i < Nthreads - 1; ++i) {
                auto block_end = block_start + block_size;
                threads[i] = std::thread(&NBodySimulation::updatePositionThread, this, block_start, block_end, dt);
                block_start = block_end;
            }
            updatePositionThread(block_start, bodies.end(), dt);
            for (auto& t : threads) t.join();

            history.push_back(bodies);
            currentTime += dt;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " ms" << std::endl;
        return history;
    }

    std::vector<std::vector<Body>> runParallelWithMutexWithTrace(double dt, size_t Nthreads) {
        std::vector<std::vector<Body>> history;
        auto startTime = std::chrono::high_resolution_clock::now();
        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return history;
        if (Nthreads > length) Nthreads = length;

        std::vector<std::mutex> force_mutexes(length);

        for (int step = 0; step < numSteps; ++step) {
            for (auto& b : bodies) b.fx = b.fy = 0.0;

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
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(computeForces, start, end);
            }
            for (auto& t : threads) t.join();

            threads.clear();
            auto updateFunc = [&](size_t start, size_t end) {
                updatePositionThread(bodies.begin() + start, bodies.begin() + end, dt);
            };

            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(updateFunc, start, end);
            }
            for (auto& t : threads) t.join();

            history.push_back(bodies);
            currentTime += dt;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation with mutex completed in " << duration << " ms" << std::endl;
        return history;
    }

    std::vector<std::vector<Body>> runParallelNoMutexWithTrace(double dt, size_t Nthreads) {
        std::vector<std::vector<Body>> history;
        auto startTime = std::chrono::high_resolution_clock::now();
        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return history;
        if (Nthreads > length) Nthreads = length;

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
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(computeForces, t, start, end);
            }
            for (auto& t : threads) t.join();

            for (size_t i = 0; i < length; ++i) {
                bodies[i].fx = bodies[i].fy = 0.0;
                for (size_t t = 0; t < Nthreads; ++t) {
                    bodies[i].fx += force_acc[i][t].first;
                    bodies[i].fy += force_acc[i][t].second;
                }
            }

            threads.clear();
            auto updateFunc = [&](size_t start, size_t end) {
                updatePositionThread(bodies.begin() + start, bodies.begin() + end, dt);
            };

            for (size_t t = 0; t < Nthreads; ++t) {
                size_t start = t * block_size;
                size_t end;
                if (t == Nthreads - 1) {
                    end = length;
                } else {
                    end = start + block_size;
                }
                threads.emplace_back(updateFunc, start, end);
            }
            for (auto& t : threads) t.join();

            history.push_back(bodies);
            currentTime += dt;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation without mutex completed in " << duration << " ms" << std::endl;
        return history;
    }

    std::vector<std::vector<Body>> runBarnesHuttWithTrace(double domainSize, double dt) {
        std::vector<std::vector<Body>> history;
        auto startTime = std::chrono::high_resolution_clock::now();
        Quadtree tree;

        for (int step = 0; step < numSteps; ++step) {
            tree.reset();
            int root = tree.createNode(0.0, 0.0, domainSize / 2.0);

            for (int i = 0; i < bodies.size(); ++i) {
                tree.insert(root, i, bodies);
            }

            for (Body& b : bodies) {
                b.fx = b.fy = 0.0;
                tree.computeForce(root, b, bodies);
            }

            tree.update_pos_vel(bodies, dt);
            history.push_back(bodies);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Barnes Hutt simulation completed in " << duration << " ms" << std::endl;
        return history;
    }


};

////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// SIMULATIONS TO TEST///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void createRandomSystem(NBodySimulation& sim, int numBodies) {
    srand(static_cast<unsigned int>(time(nullptr)));
    // central massive body (like a star) 
    sim.newBody(Body(1.0e30, 0.0, 0.0, 0.0, 0.0)); 
    
    // add random bodies
    for (int i = 1; i < numBodies; ++i) {
        double angle = (rand() % 1000) * 2.0 * M_PI / 1000;
        double distance = 1.0e11 + (rand() % 1000) * 1.0e9;
        double mass = 1.0e23 + (rand() % 1000) * 1.0e22;
        
        double x = distance * cos(angle);
        double y = distance * sin(angle);
        
        // calculate circular orbit velocity
        double v = std::sqrt(G * 1.0e30 / distance)*10000000;
        double vx = -v * sin(angle);
        double vy = v * cos(angle);
        
        sim.newBody(Body(mass, x, y, vx, vy));
    }
}

void createSolarSystem(NBodySimulation& sim, double angle_offset_radians) {
    // SUN
    sim.newBody(Body(1.9885e30, 0.0, 0.0, 0.0, 0.0)); // Sun

    struct Planet {
        double mass; // in kg
        double distance; // in meters
        const char* name; 
    };
    // sources for the orbital radius values (in meters): source: NASA, JPL, and Wikipedia
    // https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    //https://ssd.jpl.nasa.gov/planets/phys_par.html
    std::vector<Planet> planets = {
        {3.3011e23,      57909227000.0,     "Mercury"},
        {4.8675e24,     108209475000.0,     "Venus"},
        {5.97237e24,    149598262000.0,     "Earth"},
        {6.4171e23,     227943824000.0,     "Mars"},
        {1.8982e27,     778340821000.0,     "Jupiter"},
        {5.6834e26,    1426666422000.0,     "Saturn"},
        {8.6810e25,    2870658186000.0,     "Uranus"},
        {1.02413e26,   4498396441000.0,     "Neptune"}
    };

    /*for (const auto& p : planets) {
        double x = p.distance;
        double y = 0.0;

        // Orbital speed: circular orbit approximation
        double v = std::sqrt(G * 1.9885e30 / p.distance); 
        double vx = 0.0;
        double vy = v;

        sim.newBody(Body(p.mass, x, y, vx, vy));
    }*/
    for (const auto& p : planets) {
        double R = p.distance;

        // position rotated by angle_offset
        double x = R * std::cos(angle_offset_radians);
        double y = R * std::sin(angle_offset_radians);

        // orbital speed
        double v = std::sqrt(G * 1.9885e30 / R);

        // velocity that is perpendicular to radius vector so it is a 90 degree rotation
        double vx = -v * std::sin(angle_offset_radians);
        double vy =  v * std::cos(angle_offset_radians);

        sim.newBody(Body(p.mass, x, y, vx, vy));
    }

}


double computeAccuracy(const std::vector<Body>& reference, const std::vector<Body>& test) {
        if (reference.size() != test.size()) return 0.0;

        double total_diff = 0.0;
        double total_ref = 0.0;

        for (size_t i = 0; i < reference.size(); ++i) {
            double dx = reference[i].x - test[i].x;
            double dy = reference[i].y - test[i].y;
            double dvx = reference[i].vx - test[i].vx;
            double dvy = reference[i].vy - test[i].vy;

            double diff = std::sqrt(dx * dx + dy * dy + dvx * dvx + dvy * dvy);
            double ref_magnitude = std::sqrt(reference[i].x * reference[i].x + reference[i].y * reference[i].y +
                                            reference[i].vx * reference[i].vx + reference[i].vy * reference[i].vy);

            total_diff += diff;
            total_ref += ref_magnitude;
        }

        if (total_ref == 0.0) return 100.0;
        double accuracy = 100.0 * (1.0 - (total_diff / total_ref));
        return std::max(0.0, accuracy); // clamp to [0, 100]
    }

double computeAccuracyOverTime(const std::vector<std::vector<Body>>& reference_steps,
                               const std::vector<std::vector<Body>>& test_steps) {
    if (reference_steps.size() != test_steps.size()) return 0.0;

    size_t time_steps = reference_steps.size();
    size_t num_bodies = reference_steps[0].size();

    double total_accuracy = 0.0;
    size_t valid_steps = 0;

    for (size_t t = 0; t < time_steps; ++t) {
        const auto& ref = reference_steps[t];
        const auto& test = test_steps[t];

        if (ref.size() != test.size()) continue;

        double step_accuracy_sum = 0.0;
        size_t valid_bodies = 0;

        for (size_t i = 0; i < num_bodies; ++i) {
            double dx = ref[i].x - test[i].x;
            double dy = ref[i].y - test[i].y;
            double dvx = ref[i].vx - test[i].vx;
            double dvy = ref[i].vy - test[i].vy;

            double diff = std::sqrt(dx*dx + dy*dy + dvx*dvx + dvy*dvy);
            double ref_mag = std::sqrt(ref[i].x*ref[i].x + ref[i].y*ref[i].y +
                                       ref[i].vx*ref[i].vx + ref[i].vy*ref[i].vy);

            if (ref_mag > 0.0) {
                double body_accuracy = 100.0 * (1.0 - diff / ref_mag);
                step_accuracy_sum += std::max(0.0, body_accuracy);
                ++valid_bodies;
            }
        }

        if (valid_bodies > 0) {
            total_accuracy += step_accuracy_sum / valid_bodies;
            ++valid_steps;
        }
    }

    if (valid_steps == 0) return 0.0;
    return total_accuracy / valid_steps;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// MAIN FUNCTION ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    /*std::ofstream timingFile("timing_results.csv", std::ios::app);
    if (!timingFile.is_open()) {
        std::cerr << "Could not open timing_results.csv" << std::endl;
        return 1;
    }

    std::ofstream accuracyFile("accuracy_results.csv", std::ios::app);
    if (!accuracyFile.is_open()) {
        std::cerr << "Could not open accuracy_results.csv" << std::endl;
        return 1;
    }*/

    int num_threads = 6; // default
    /*if (argc > 1) {
        num_threads = std::stoi(argv[1]);
    }*/
    int num_bodies = 50;
    if (argc > 1) {
        num_bodies = std::stoi(argv[1]);
    }

    // simulation with time step and total time
    NBodySimulation simulation_sequential(1, 50); // 1 second time step, 20 total
    //NBodySimulation simulation_sequential(3600 * 24, 3600 * 24 * 365); // 1 day timestep, simulate 1Y

    // system of bodies for simulation
    createRandomSystem(simulation_sequential, num_bodies); 
    //createSolarSystem(simulation_sequential, 0.0);
    // createSolarSystem(simulation_sequential, M_PI/4.0); // centered solar system with planets circlinng the sun (45 degree angle)

    std::vector<Body> initial_bodies = simulation_sequential.getBodies();
    //std::cout << "TOTAL NUMBER OF BODIES" << initial_bodies.size() << "\n";
    
    // Run sequential simulation
    //long seq_time = simulation_sequential.runSequential();
    std::vector<std::vector<Body>> trace_seq = simulation_sequential.runSequentialWithTrace();
    std::vector<Body> sequential_result = simulation_sequential.getBodies();

    NBodySimulation simulation_parallel(1, 50);
    //NBodySimulation simulation_parallel(3600 * 24, 3600 * 24 * 365); // 1 day timestep, simulate 1Y
    simulation_parallel.setBodies(initial_bodies);
    //long par_time = simulation_parallel.runParallel(1.0, num_threads);
    std::vector<std::vector<Body>> trace_par = simulation_parallel.runParallelWithTrace(1.0, num_threads);
    std::vector<Body> parallel_result = simulation_parallel.getBodies();


    NBodySimulation simulation_parallel_mutex(1, 50);
    //NBodySimulation simulation_parallel_mutex(3600 * 24, 3600 * 24 * 365); // 1 day timestep, simulate 1Y
    simulation_parallel_mutex.setBodies(initial_bodies);
    //long mutex_time = simulation_parallel_mutex.runParallelWithMutex(1.0, num_threads);
    std::vector<std::vector<Body>> trace_mutex = simulation_parallel_mutex.runParallelWithMutexWithTrace(1.0, num_threads);
    std::vector<Body> mutex_result = simulation_parallel_mutex.getBodies();


    NBodySimulation simulation_parallel_nomutex(1, 50);
    //NBodySimulation simulation_parallel_nomutex(3600 * 24, 3600 * 24 * 365); // 1 day timestep, simulate 1Y
    simulation_parallel_nomutex.setBodies(initial_bodies);
    //long nomutex_time = simulation_parallel_nomutex.runParallelNoMutex(1.0, num_threads);
    std::vector<std::vector<Body>> trace_nomutex = simulation_parallel_nomutex.runParallelNoMutexWithTrace(1.0, num_threads);
    std::vector<Body> nomutex_result = simulation_parallel_nomutex.getBodies();

    //NBodySimulation for Barnes Hutt 
    NBodySimulation simulation_barneshutt(1, 50);
    //NBodySimulation simulation_barneshutt(3600 * 24, 3600 * 24 * 365);
    simulation_barneshutt.setBodies(initial_bodies);
    double domainSize = 1e13;// should be large enough to include all bodies
    //simulation_barneshutt.runBarnesHutt(domainSize, 1.0); // timestep = 1.0
    std::vector<std::vector<Body>> trace_barneshutt = simulation_barneshutt.runBarnesHuttWithTrace(domainSize, 1.0);
    std::vector<Body> barneshutt_result = simulation_barneshutt.getBodies();

    //timingFile << num_threads << "," << seq_time << "," << par_time << "," << mutex_time << "," << nomutex_time << "\n";
    //timingFile.close();

    /*bool same_results = true;

    if (sequential_result.size() != parallel_result.size()) {
        same_results = false;
    } 
    
    else {
        for (size_t i=0; i < sequential_result.size(); ++i) {
            double dx = std::abs(sequential_result[i].x - parallel_result[i].x);
            double dy = std::abs(sequential_result[i].y - parallel_result[i].y);
            double dvx = std::abs(sequential_result[i].vx - parallel_result[i].vx);
            double dvy = std::abs(sequential_result[i].vy - parallel_result[i].vy);
            std::cout << "dx sequential" << dx << "\n";
            std::cout << "dy seq " << dy << "\n";
            std::cout << "dvx seq" << dvx << "\n";
            std::cout << "dvy seq" << dvy << "\n";

            //if (dx > 1e-6 || dy > 1e-6 || dvx > 1e-6 || dvy > 1e-6) {
            if (dx > 1e9 || dy > 1e8 || dvx > 1e1 || dvy > 1e1) {
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
    */


    // Test if both simulations grant the same result (ChatGPT helped me debug the code by adding the comparison of sizes before checking values and added the 'break')
    auto compareResults = [](const std::vector<Body>& a, const std::vector<Body>& b) {
        if (a.size() != b.size()) return false;

        for (size_t i = 0; i < a.size(); ++i) {
            double dx = std::abs(a[i].x - b[i].x);
            double dy = std::abs(a[i].y - b[i].y);
            double dvx = std::abs(a[i].vx - b[i].vx);
            double dvy = std::abs(a[i].vy - b[i].vy);
            //std::cout << "dx" << dx << "\n";
            //std::cout << "dy" << dy << "\n";
            //std::cout << "dvx" << dvx << "\n";
            //std::cout << "dvy" << dvy << "\n";


            if (dx > 1e-1 || dy > 1e-1 || dvx > 1e-1 || dvy > 1e-1) {
            //if (dx > 1e9 || dy > 1e8 || dvx > 1e1 || dvy > 1e1) {
                return false;
            }
        }
        return true;
    };
    if (compareResults(sequential_result, parallel_result)) {
       std::cout << "OK." << std::endl;
    } else {
        std::cout << "Parallel and sequential simulations grant different results" << std::endl;
    }
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
    if (compareResults(sequential_result, barneshutt_result)) {
        std::cout << "OK." << std::endl;
    } else {
        std::cout << "Sequential and Barnes-Hutt simulations grant different results" << std::endl;
    }
    
    

    double acc_parallel = computeAccuracyOverTime(trace_seq, trace_par);
    double acc_mutex = computeAccuracyOverTime(trace_seq, trace_mutex);
    double acc_nomutex = computeAccuracyOverTime(trace_seq, trace_nomutex);
    double acc_barneshutt = computeAccuracyOverTime(trace_seq, trace_barneshutt);

    // Output to CSV
    std::ofstream out("accuracy_results.csv", std::ios::app);
    out << num_bodies << "," << acc_parallel << "," << acc_mutex << "," << acc_nomutex << "," << acc_barneshutt << "\n";
    out.close();

    /*std::ofstream accFile("accuracy_results.csv", std::ios::app);
    if (accFile.is_open()) {
        accFile << initial_bodies.size() << "," 
                << acc_parallel << "," 
                << acc_mutex << "," 
                << acc_nomutex << "," 
                << acc_barneshutt << "\n";
        accFile.close();
    }*/
    //accuracyFile << initial_bodies.size() << "," << acc_parallel << "," << acc_mutex << "," << acc_nomutex << "," << acc_barneshutt << "\n";
    //accuracyFile.close();

    ////// Compare results of sequential and barnes hutt

    /*same_results = true;

    if (sequential_result.size() != barneshutt_result.size()) {
        same_results = false;
    } 
    
    else {
        for (size_t i=0; i < sequential_result.size(); ++i) {
            double dx = std::abs(sequential_result[i].x - barneshutt_result[i].x);
            double dy = std::abs(sequential_result[i].y - barneshutt_result[i].y);
            double dvx = std::abs(sequential_result[i].vx - barneshutt_result[i].vx);
            double dvy = std::abs(sequential_result[i].vy - barneshutt_result[i].vy);
            std::cout << "dx sequential" << dx << "\n";
            std::cout << "dy seq " << dy << "\n";
            std::cout << "dvx seq" << dvx << "\n";
            std::cout << "dvy seq" << dvy << "\n";

            if (dx > 1e-6 || dy > 1e-6 || dvx > 1e-6 || dvy > 1e-6) {
                same_results = false;
                break;
            }
        }
    }

    if (same_results) {
        std::cout << "OK." << std::endl;
    } else {
        std::cout << "Sequential and Barnes-Hutt simulations grant different results" << std::endl;
    }
    return 0;*/
}

/*
When compiling:  make ./simple_approach_test
When running: ./simple_approach_test

Run the Python file in venv environment:
python make_frames.py

Back in terminal:
magick convert -delay 20 -loop 0 frames1/frame_*.png simulation1.gif


1) When we only had runSequential
Sequential simulation completed in 7 ms

2) When we added parallelization
Sequential simulation completed in 7 ms
Parallel simulation completed in 13 ms

3) When we avoided computing Fij twice
Sequential simulation completed in 9 ms
Parallel simulation completed in 7 ms

*/
