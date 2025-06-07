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
#include <random>
#include <complex> 
#include <fstream>
#include <algorithm>
#include <future>
#include <iomanip>

#define G 6.67430e-11 // Gravity constant
#define NOT_ZERO 1e-9 // constant to avoid division by zero
#define THETA 0.1 // Barnes hut approximation parameter

// constants for the FMM method :
#define P 20 // number of multipole/local expansion terms
#define MAX_BODIES_PER_CELL 2. // maximum number of bodies per leaf cell

double sqr(double x) {
    return x*x;
}

using Complex = std::complex<double>;

// _______________________________________________________________________________________________________
// Structure to represent a body in 2D space 
class Body {
public:
    double mass;     // of the body
    double x, y;     // position coordinates
    double vx, vy;   // velocity components
    double fx, fy;   // Force components
    Complex position;


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

    std::complex<double> get_position() const {
        std::complex<double> position(x, y);
        return position;}
    
    std::complex<double> get_velocity() const {
        std::complex<double> velocity(vx, vy);
        return velocity;}
    
    std::complex<double> get_force() const {
        std::complex<double> force(fx, fy);
        return force;}
    
};

// _______________________________________________________________________________________________________

class BBox{ // class bounding box for the FMM method
public: 
    double xmin; 
    double xmax;
    double ymin; 
    double ymax;

    BBox(double xmin = 0, double xmax = 1, double ymin = 0, double ymax = 1)
        : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax) {}
    
    std::complex<double> get_center() const {
        std::complex<double> center((xmin + xmax) / 2, (ymin + ymax) / 2);
         return center;}

    std::complex<double> get_position(Body& body) const {
        std::complex<double> position(body.x, body.y);
        return position;}

    double size() const { return std::max(xmax - xmin, ymax - ymin); }
    
    bool contains(const std::complex<double> & point) const {
        return point.real() >= xmin && point.real() <= xmax &&
               point.imag() >= ymin && point.imag() <= ymax;
    }
};

class FMMNode {
public:
    BBox bbox;
    std::vector<int> body_id;
    std::vector<std::unique_ptr<FMMNode>> children;
    bool is_leaf;
    
    // FMM coefficients
    std::vector<std::complex<double>> multipole_coeffs;  // a_k coefficients
    std::vector<std::complex<double>> local_coeffs;      // b_k coefficients
    std::complex<double> center;
    
    // Parallel processing support
    mutable std::mutex node_mutex;
    
    FMMNode(const BBox& box) : bbox(box), is_leaf(true), center(box.get_center()) {
        multipole_coeffs.resize(P, 0.0);
        local_coeffs.resize(P, 0.0);
    }
    
    void subdivide() {
        if (is_leaf && body_id.size() > MAX_BODIES_PER_CELL) {
            is_leaf = false;
            children.resize(4);
            
            double midx = (bbox.xmin + bbox.xmax) / 2;
            double midy = (bbox.ymin + bbox.ymax) / 2;
            
            // divide the previous box by 2 
            children[0] = std::make_unique<FMMNode>(BBox(bbox.xmin, midx, bbox.ymin, midy)); // SW
            children[1] = std::make_unique<FMMNode>(BBox(midx, bbox.xmax, bbox.ymin, midy)); // SE
            children[2] = std::make_unique<FMMNode>(BBox(bbox.xmin, midx, midy, bbox.ymax)); // NW
            children[3] = std::make_unique<FMMNode>(BBox(midx, bbox.xmax, midy, bbox.ymax)); // NE
        }
    }
};

// Parallel FMM solver
class ParallelFMMSolver {
private:
    std::vector<Body> bodies;
    std::unique_ptr<FMMNode> root;
    BBox box;
    unsigned int Nthreads;


public:
    ParallelFMMSolver(const std::vector<Body>& bodies, const BBox& box, unsigned int Nthreads = std::thread::hardware_concurrency()) 
        : bodies(bodies), box(box), Nthreads(Nthreads) {
        root = std::make_unique<FMMNode>(box);
    }

    void setBodies(const std::vector<Body>& newBodies) { 
        bodies = newBodies; 
    }

    const std::vector<Body>& getBodies() { 
        return bodies; 
    }
    
    void buildTree() {
        root = std::make_unique<FMMNode>(box);
        
        // add all bodies to root
        for (int i = 0; i < bodies.size(); ++i) {
            root->body_id.push_back(i);
        }
        
        buildTreeRecursive(root.get());
    }
    
    void buildTreeRecursive(FMMNode* node) {
        if (node->body_id.size() <= MAX_BODIES_PER_CELL) {
            return; // leaf node
        }
        
        node->subdivide();
        
        // distribute bodies to children
        for (int idx : node->body_id) {
            for (auto& child : node->children) {
                std::complex<double> pos(bodies[idx].x, bodies[idx].y);
                if (child->bbox.contains(pos)) {
                    child->body_id.push_back(idx);
                    break;
                }
            }
        }
        
        node->body_id.clear(); 
        
        // Recursively build children
        for (auto& child : node->children) {
            if (!child->body_id.empty()) {
                buildTreeRecursive(child.get());
            }
        }
    }
    void computeMultipoleParaAUX(FMMNode* node) {
        if (node->is_leaf) {
            node->multipole_coeffs[0] = 0.0;
            for (int idx : node->body_id) {
                const Body& p = bodies[idx];
                Complex z = p.get_position() - node->center;

                node->multipole_coeffs[0] += p.mass; 

                Complex z_pow = 1.0;
                for (int k = 1; k < P; ++k) {
                    z_pow *= -z;
                    node->multipole_coeffs[k] += p.mass * z_pow / double(k);
                }
            }
        } else {
            std::fill(node->multipole_coeffs.begin(), node->multipole_coeffs.end(), 0.0);

            for (auto& child : node->children) {
                if (child && !child->body_id.empty()) {
                    Complex z0 = child->center - node->center;
                    translateMultipole(child->multipole_coeffs, node->multipole_coeffs, z0);
                }
            }
        }
    }

    void computeMultipoleExpansionsParallel() {
        std::vector<std::thread> threads;
        std::vector<FMMNode*> nodesToProcess;

        // collect all nodes that need processing (don't take empty nodes)
        collectNodesForMultipoleExpansion(root.get(), nodesToProcess);

        // compute multipole expansions in parallel
        for (auto* node : nodesToProcess) {
            threads.emplace_back([this, node]() {
                computeMultipoleParaAUX(node);
            });
        }
        for (auto& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
    
    void collectNodesForMultipoleExpansion(FMMNode* node, std::vector<FMMNode*>& nodes) {
        if (node->is_leaf || !node->body_id.empty()) {
            nodes.push_back(node);
        } else {
            for (auto& child : node->children) {
                if (child) {
                    collectNodesForMultipoleExpansion(child.get(), nodes);
                }
            }
        }
    }
    
    void translateMultipole(const std::vector<Complex>& source, std::vector<Complex>& target, Complex z0) {
        target[0] += source[0];
        
        for (int k = 1; k < P; ++k) {
            Complex sum = source[k];
            Complex z0_pow = 1.0;
            
            for (int j = 1; j < k; ++j) {
                z0_pow *= z0;
                sum += source[k-j] * z0_pow / double(j);
            }
            
            target[k] += sum;
        }
    }
    
    void computeForces() { 
        // reset forces before computing the forces to avoid accumulated forces over time
        for (auto& p : bodies) {
            p.fx = 0.0;
            p.fy = 0.0;
        }
        computeForcesRecursive(root.get(), root.get());
    }
    
    void computeForcesRecursive(FMMNode* source, FMMNode* target) {
        if (!source || !target || source->body_id.empty() || target->body_id.empty()) {
            return;
        }
        
        Complex r = target->center - source->center;
        double distance = abs(r);
        double source_size = source->bbox.size();

        if (source != target && source_size / distance < THETA) {
            // Use multipole approximation
            applyMultipoleForce(source, target);
        } else if (source->is_leaf && target->is_leaf) {
            // direct bodies interaction
            directInteraction(source, target);
        } else {
            // recurse to children
            if (!source->is_leaf) {
                for (auto& child : source->children) {
                    computeForcesRecursive(child.get(), target);
                }
            } else {
                for (auto& child : target->children) {
                    computeForcesRecursive(source, child.get());
                }
            }
        }
    }
    
    void applyMultipoleForce(FMMNode* source, FMMNode* target) {
        Complex z0 = target->center - source->center;
        
        for (int target_idx : target->body_id) {
            Body& p = bodies[target_idx];
            Complex z = p.get_position() - source->center;
            double r = std::abs(z);
            if (r < 1e-10) continue;

            double force_mag = G * p.mass * source->multipole_coeffs[0].real() / (r * r);
            Complex force = z / r; 

            p.fx += G * p.mass * force.real();
            p.fy += G * p.mass * force.imag();
        }
    }
    
    void directInteraction(FMMNode* source, FMMNode* target) {
        for (int i : source->body_id) {
            for (int j : target->body_id) {
                if (i != j) {
                    const Body& pi = bodies[i];
                    Body& pj = bodies[j];
                    
                    Complex r = pj.get_position() - pi.get_position();
                    double dist = abs(r);
                    
                    if (dist > 1e-10) { // avoid division by zero
                        Complex force_dir = r / dist;
                        double force_mag = G * pi.mass * pj.mass / (dist * dist);
                        pj.fx += force_mag * force_dir.real();
                        pj.fy += force_mag * force_dir.imag();
                    }
                }
            }
        }
    }
    
    void updateParticles(double dt) {
        for (auto& p : bodies) {
            p.updateVelocity(dt);
            p.updatePosition(dt);
        }
    }
    
    void simulate(double dt, int steps) {
        auto start = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < steps; ++step) {
            buildTree();
            computeMultipoleExpansionsParallel();
            computeForces();
            updateParticles(dt);
            
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::ofstream filebis("results/results_parallel_fmm.csv", std::ios_base::app);
        filebis << bodies.size() << ',' << duration << ',';
        filebis << '\n'; 
        filebis.close();
        
        std::cout << "Parallel FMM simulation completed in " << duration << " ms" << std::endl;
    }

    const std::vector<Body>& getParticles() const { return bodies; }
};


class FMMSolver {
private:
    std::vector<Body> bodies;
    std::unique_ptr<FMMNode> root;
    BBox box;
    
    
public:
    FMMSolver(const std::vector<Body>& bodies, const BBox& box) 
        : bodies(bodies), box(box) {
        root = std::make_unique<FMMNode>(box);
    }

    void setBodies(const std::vector<Body>& newBodies) { 
        bodies = newBodies; 
    }

    const std::vector<Body>& getBodies() { 
        return bodies; 
    }
    
    void buildTree() {
        root = std::make_unique<FMMNode>(box);
        
        // add all bodies to root
        for (int i = 0; i < bodies.size(); ++i) {
            root->body_id.push_back(i);
        }
        
        buildTreeRecursive(root.get());
    }
    
    void buildTreeRecursive(FMMNode* node) {
        if (node->body_id.size() <= MAX_BODIES_PER_CELL) {
            return; // leaf node
        }
        
        node->subdivide();
        
        // distribute bodies to children
        for (int idx : node->body_id) {
            for (auto& child : node->children) {
                std::complex<double> pos(bodies[idx].x, bodies[idx].y);
                if (child->bbox.contains(pos)) {
                    child->body_id.push_back(idx);
                    break;
                }
            }
        }
        
        node->body_id.clear(); 
        
        // Recursively build children
        for (auto& child : node->children) {
            if (!child->body_id.empty()) {
                buildTreeRecursive(child.get());
            }
        }
    }
    
    void computeMultipoleExpansions() {
        computeMultipoleRecursive(root.get());
    }
    
    void computeMultipoleRecursive(FMMNode* node) {
        if (node->is_leaf) {
            // compute multipole expansion for leaf node
            node->multipole_coeffs[0] = 0.0;
            for (int idx : node->body_id) {
                const Body& p = bodies[idx];
                Complex z = p.get_position() - node->center;
                
                node->multipole_coeffs[0] += p.mass; // total mass
                
                Complex z_pow = 1.0;
                for (int k = 1; k < P; ++k) {
                    z_pow *= -z; 
                    node->multipole_coeffs[k] += p.mass * z_pow / double(k);
                }
            }
        } else {
            // compute multipole expansion from children
            fill(node->multipole_coeffs.begin(), node->multipole_coeffs.end(), 0.0);
            
            for (auto& child : node->children) {
                if (child && !child->body_id.empty()) {
                    computeMultipoleRecursive(child.get());
                    
                    // translate child's multipole expansion to parent center
                    Complex z0 = child->center - node->center;
                    translateMultipole(child->multipole_coeffs, node->multipole_coeffs, z0);
                }
            }
        }
    }
    
    void translateMultipole(const std::vector<Complex>& source, std::vector<Complex>& target, Complex z0) {
        target[0] += source[0];
        
        for (int k = 1; k < P; ++k) {
            Complex sum = source[k];
            Complex z0_pow = 1.0;
            
            for (int j = 1; j < k; ++j) {
                z0_pow *= z0;
                sum += source[k-j] * z0_pow / double(j);
            }
            
            target[k] += sum;
        }
    }
    
    void computeForces() {
        //Reset forces
        for (auto& p : bodies) {
            p.fx = 0.0;
            p.fy = 0.0;
        }
        
        computeForcesRecursive(root.get(), root.get());
    }
    
    void computeForcesRecursive(FMMNode* source, FMMNode* target) {
        if (!source || !target || source->body_id.empty() || target->body_id.empty()) {
            return;
        }
        
        Complex r = target->center - source->center;
        double distance = abs(r);
        double source_size = source->bbox.size();

        if (source != target && source_size / distance < THETA) {
            // use the multipole approximation
            applyMultipoleForce(source, target);
        } else if (source->is_leaf && target->is_leaf) {
            // Direct bodies interaction
            directInteraction(source, target);
        } else {
            // recurse to children
            if (!source->is_leaf) {
                for (auto& child : source->children) {
                    computeForcesRecursive(child.get(), target);
                }
            } else {
                for (auto& child : target->children) {
                    computeForcesRecursive(source, child.get());
                }
            }
        }
    }
    
    void applyMultipoleForce(FMMNode* source, FMMNode* target) {
        Complex z0 = target->center - source->center;
        
        for (int target_idx : target->body_id) {
            Body& p = bodies[target_idx];
            Complex z = p.get_position() - source->center;
            double r = std::abs(z);
            if (r < 1e-10) continue;

            double force_mag = G * p.mass * source->multipole_coeffs[0].real() / (r * r);
            Complex force = z / r; 

            p.fx += G * p.mass * force.real();
            p.fy += G * p.mass * force.imag();
        }
    }
    
    void directInteraction(FMMNode* source, FMMNode* target) {
        for (int i : source->body_id) {
            for (int j : target->body_id) {
                if (i != j) {
                    const Body& pi = bodies[i];
                    Body& pj = bodies[j];
                    
                    Complex r = pj.get_position() - pi.get_position();
                    double dist = abs(r);
                    
                    if (dist > 1e-10) { // Avoid division by zero
                        Complex force_dir = r / dist;
                        double force_mag = G * pi.mass * pj.mass / (dist * dist);
                        pj.fx += force_mag * force_dir.real();
                        pj.fy += force_mag * force_dir.imag();
                    }
                }
            }
        }
    }
    
    void updateParticles(double dt) {
        for (auto& p : bodies) {
            // Leapfrog integration
            p.updateVelocity(dt);
            p.updatePosition(dt);
        }
    }
    
    void simulate(double dt, int steps) {
        auto start = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < steps; ++step) {
            buildTree();
            computeMultipoleExpansions();
            computeForces();
            updateParticles(dt);
            
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::ofstream filebis("results/results_sequential_fmm.csv", std::ios_base::app);
        filebis << bodies.size() << ',' << duration << ',';
        filebis << '\n'; 
        filebis.close();
        
        std::cout << "Sequential FMM simulation completed in " << duration << " ms" << std::endl;
    }

    const std::vector<Body>& getParticles() const { return bodies; }
};

// _______________________________________________________________________________________________________

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
        std::ofstream filebis("results/results_sequential.csv", std::ios_base::app);
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
        
        filebis << bodies.size() << ','<< duration <<',';
        filebis << '\n'; 
        filebis.close();
        
        std::cout << "Sequential simulation completed in " << duration << " milliseconds" << std::endl;
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

        std::ofstream filebis("results/results_parallel.csv", std::ios_base::app);
        // End timer
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " ms" << std::endl;
        
        filebis << Nthreads << ',' << bodies.size() << ','<< duration <<',';
        filebis << '\n'; 
        filebis.close();
        

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
        std::ofstream filebis("results/results_withmutex.csv", std::ios_base::app);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " ms" << std::endl;

        filebis << Nthreads << ','<< bodies.size() << ','<< duration <<',';
        filebis << '\n'; 
        filebis.close();
    }


    void runParallelNoMutex(double dt, size_t Nthreads) {
        auto startTime = std::chrono::high_resolution_clock::now();

        size_t length = bodies.size();
        if (length == 0 || Nthreads == 0) return;
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

            currentTime += dt;
        }
        std::ofstream filebis("results/results_nomutex.csv", std::ios_base::app);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Parallel simulation completed in " << duration << " microseconds" << std::endl;

        filebis << Nthreads << ',' << bodies.size() << ','<< duration <<',';
        filebis << '\n'; 
        filebis.close();
    }
};

// _______________________________________________________________________________________________________

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


// _______________________________________________________________________________________________________
int main() {
    // Set random seed
    srand(static_cast<unsigned int>(time(nullptr)));
    
    double num_bodies = 2000; 
    double timestep = 1;
    double total_time = 20;
    int num_threads = std::thread::hardware_concurrency();
    std::cout << "num_threads: " << num_threads << std::endl;

    for (int i =0 ; i< 1; i++){
        
        std::cout << "N-Body simulation"<< std::endl;
        std::cout << "number of bodies: " << num_bodies << std::endl;

        // Create simulation with time step and total time
        NBodySimulation simulation_sequential(timestep, total_time); // 1 second time step, 20 total

        // Create a system of bodies 
        createRandomSystem(simulation_sequential, num_bodies); 
        std::vector<Body> initial_bodies = simulation_sequential.getBodies();
        

        // Run sequential simulation
        simulation_sequential.runSequential();
        std::vector<Body> sequential_result = simulation_sequential.getBodies();


        // run paralllel
        NBodySimulation simulation_parallel(1, total_time);
        simulation_parallel.setBodies(initial_bodies);
        simulation_parallel.runParallel(1.0, num_threads);
        std::vector<Body> parallel_result = simulation_parallel.getBodies();

        // run parallel with mutex 
        NBodySimulation simulation_parallel_mutex(1, total_time);
        simulation_parallel_mutex.setBodies(initial_bodies);
        simulation_parallel_mutex.runParallelWithMutex(1.0, num_threads);
        std::vector<Body> mutex_result = simulation_parallel_mutex.getBodies();


        // run parallel without mutex
        NBodySimulation simulation_parallel_nomutex(1, total_time);
        simulation_parallel_nomutex.setBodies(initial_bodies);
        simulation_parallel_nomutex.runParallelNoMutex(1.0, num_threads);
        std::vector<Body> nomutex_result = simulation_parallel_nomutex.getBodies();

        auto compareResults = [](const std::vector<Body>& a, const std::vector<Body>& b, int c, double num_bodies) {

            std::ofstream file( "results/" + std::to_string(c) + "_accuracy.csv",  std::ios_base::app);
            if (a.size() != b.size()) return false;
            bool result = true;
            for (size_t i = 0; i < a.size(); ++i) {
                double dx = std::abs(a[i].x - b[i].x);
                //std::cout << "dx : " << dx << std::endl;
                double dy = std::abs(a[i].y - b[i].y);
                //std::cout << "dy : " << dy << std::endl;
                double dvx = std::abs(a[i].vx - b[i].vx);
                //std::cout << "dvx : " << dvx << std::endl;
                double dvy = std::abs(a[i].vy - b[i].vy);
                //std::cout << "dvy : " << dvy << std::endl;
                //std::cout << "\n"<< std::endl;
                file << "num_bodies_"<< num_bodies << ', '<< "bodyi_" << i <<  ',' << std::fixed << std::setprecision(5) << dx << ',' << std::fixed << std::setprecision(5) << dy << ',' << std::fixed << std::setprecision(5) << dvx << ',' << std::fixed << std::setprecision(5) << dvy;
                file <<  '\n';
                
                if (dx > 1e-1 || dy > 1e-1 || dvx > 1e-1 || dvy > 1e-1) {
                    result = false;
                }                
            }
            file.close();
            return result;
        };

        //Test if both simulations grant the same result (ChatGPT helped me debug the code by adding the comparison of sizes before checking values and added the 'break')
        std::cout << "compare results of parallel and sequential "<< std::endl;
        std::cout << "\n";
        ////////////// 
        if (compareResults(sequential_result, parallel_result,0, num_bodies)) {
            std::cout << "OK." << std::endl;
            std::cout << "\n";
        } else {
            std::cout << "Parallel and sequential simulations grant different results" << std::endl;
            std::cout << "\n";
        }

        // Compare results of parallel and mutex 
        std::cout << "compare results of sequential and mutex "<< std::endl;
        if (compareResults(sequential_result, mutex_result, 1, num_bodies)) {
            std::cout << "OK." << std::endl;
            std::cout << "\n";
        } else {
            std::cout << "Mutex and sequential simulations grant different results" << std::endl;
            std::cout << "\n";
        }

        std::cout << "compare results of sequential and no mutex "<< std::endl;
        if (compareResults(sequential_result, nomutex_result, 2, num_bodies)) {
            std::cout << "OK." << std::endl;
            std::cout << "\n";
        } else {
            std::cout << "No mutex version and sequential simulations grant different results" << std::endl;
            std::cout << "\n";
        }


        /* Run FMM method */
        std::cout << "testing Finite Multiple Method FMM, L. Greengard and V. Rokhlin"<< std::endl;


        double max_coord = 0;
        for (const auto& body : initial_bodies) {
            max_coord = std::max(max_coord, std::max(std::abs(body.x), std::abs(body.y)));
        }
        max_coord *= 1.2; // Add some margin
        
        BBox box(-max_coord, max_coord, -max_coord, max_coord);
        
        // FIXED: Use the same initial_bodies for FMM
        FMMSolver solver(initial_bodies, box);
        solver.setBodies(initial_bodies);
        
        // Run FMM simulation with same parameters
        solver.simulate(timestep, total_time / timestep);
        std::vector<Body> results_fmm = solver.getParticles();


        // Compare results
        std::cout << "\nCompare results FMM and sequential" << std::endl;
        std::cout << "Sequential results size: " << sequential_result.size() << std::endl;
        std::cout << "FMM results size: " << results_fmm.size() << std::endl;


        std::cout << "compare results FMM and sequential" << std::endl;
        if (compareResults(sequential_result, results_fmm, 3, num_bodies)) {
            std::cout << "OK." << std::endl;
            std::cout << "\n";
        } else {
            std::cout << "FMM method and sequential simulations grant different results" << std::endl;
            std::cout << "\n";
        }
        
        // test parallel FMM
        ParallelFMMSolver solverparallel(initial_bodies, box, num_threads);
        solverparallel.setBodies(initial_bodies);
        
        // Run FMM simulation with same parameters
        solverparallel.simulate(timestep, total_time / timestep);
        std::vector<Body> results_fmm_para = solverparallel.getParticles();


        // Compare results
        std::cout << "\nCompare results FMM Parallel and sequential" << std::endl;
        std::cout << "Sequential results size: " << sequential_result.size() << std::endl;
        std::cout << "FMM parallel results size: " << results_fmm_para.size() << std::endl;


        std::cout << "compare results FMM parallel and sequential" << std::endl;
        if (compareResults(sequential_result, results_fmm_para, 4, num_bodies)) {
            std::cout << "OK." << std::endl;
            std::cout << "\n";
        } else {
            std::cout << "FMM Parallel method and sequential simulations grant different results" << std::endl;
            std::cout << "\n";
        }

        num_bodies *=10;
    }
    return 0;
}

/*
When compiling:  g++ -std=c++11 parallel_body_fmmm.cpp -o bodies
When running: ./bodies


*/
