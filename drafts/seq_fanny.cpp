/* Projet 1 : N-body simulations */
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

// Universal gravitational constant (G) in m³/(kg·s²)
const double G = 6.67430e-11;

double square(double i){
    return i*i;
}

std::vector<double> forces_brute_force(const std::vector<double>& masses,const  std::vector<std::pair<double, double>>& coordinates, const std::vector<std::pair<double, double>>& velocities, size_t N ) {
    std::vector<double> forces(N, 0.0) ;
    for (size_t i=0; i<N; i++) {
        for (size_t j=0; j<N; j++) {
            if (i==j){continue;}

            double dx = coordinates[i].first-coordinates[j].first;
            double dy = coordinates[i].second-coordinates[j].second;
            if (square(dx) + square(dy)==0){continue;}
            double f= G * (masses[i]* masses[j]) / (square(dx) + square(dy));
            forces[i] += f;
        }
    }
    return forces;
}

int main() {

    const int N = 10;  
    std::vector<double> masses(N);
    std::vector<std::pair<double, double>> coordinates(N);
    std::vector<std::pair<double, double>> velocities(N);

    // Example: Initialize the vectors with some values
    // Here we simply fill them with sample data for demonstration.
    for (int i = 0; i < N; ++i) {
        masses[i] = 50.0 + i;                 // For example, m = 1.0, 2.0, 3.0, ..., etc.
        coordinates[i] = {i * 1.0, i * 2.0}; // (xi, yi) = (0, 0), (1, 2), (2, 4), ...
        velocities[i] = {i * 0.1, i * 0.2};    // (ui, vi) = (0, 0), (0.1, 0.2), (0.2, 0.4), ...
    }

    std::vector<double> forces = forces_brute_force(masses, coordinates, velocities, N);


    // Print the values to verify initialization
    for (int i = 0; i < N; ++i) {
        std::cout << "Particle " << i << ":\n";
        std::cout << "  Mass = " << masses[i] << "\n";
        std::cout << "  Coordinates = (" << coordinates[i].first << ", " << coordinates[i].second << ")\n";
        std::cout << "  Velocity = (" << velocities[i].first << ", " << velocities[i].second << ")\n";
        std::cout << "--------------------------\n";
        std::cout << "Force on object " << i << ": " << forces[i] << " N" << std::endl;
    }

    return 0;
}
