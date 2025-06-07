#include <iostream>
#include "main.cpp"
#include <chrono>


/* This code tests and analyses the different times and output of the sequential and parallelise algorithms */

int main(){

    std::cout << "This code analyse the different times and output of the sequential and parallelise algorithms "<<std::endl;

    std::cout << "Testing correctness" << std::endl;
    const size_t NUM_TESTS = 1000;
    size_t tests_passed = 0;

    // N-bodies definitions 
    // format : Body(double mass, double x, double y, double vX, double vY)
    Body earth = Body(5.97237e24,0, 0, 0, 0);
    earth.idBody = 0;
    Body moon = Body(7.342e22, 384400000.0, 0, 0, 0);
    moon.idBody = 1;
    
    std::vector<Body> bodies;
    bodies.push_back(earth);
    bodies.push_back(moon);
    size_t N = bodies.size();

    // sequential 
    
    std::vector<std::vector<double>> results(N, std::vector<double>(N, 0));
    auto startTime = std::chrono::high_resolution_clock::now();

    results = fixed_seq_forces(bodies, N);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();

    std::cout << "Force on object earth and moon : " << results[0][1] << " N" << std::endl;
    std::cout << "Sequential simulation completed in " << duration << " ms" << std::endl;



    return 0;
}
