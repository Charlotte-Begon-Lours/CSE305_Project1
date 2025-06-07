#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

/* Define the constants */
const double G = 6.67430e-11;
const double NOT_ZERO = 1e-9;


double square(double x){
    return x*x;
}

class Body{
public:
    int idBody;
    double mass;
    double x, y;
    double vx, vy;
    double aX, aY;
    double fx, fy;
    Body(double mass, double x, double y, double vx, double vy)
    : mass(mass), x(x), y(y), vx(vx), vy(vy), aX(0), aY(0){}

    void update_velocity(double dt){
        vx += dt * aX;
        vy += dt * aY;
        // vx += dt * fx / mass;
        // vy += dt * fy / mass;
    }
    void update_position(double dt){
        x += dt * vx;
        y += dt * vy;
    }

    void Force_update(const Body& other) { // force of another body on this body
        double dx = other.x - x;
        double dy = other.y - y;
        double dist = std::sqrt(square(dx) + square(dy));
        
        // To avoid division by zero
        dist = std::max(dist, NOT_ZERO);

        // Newton's law:
        double force = G * mass * other.mass / (dist * dist);
        
        // add these forces to the total exerted force
        fx += force * dx / dist;
        fy += force * dy / dist;
    }

};

// store the forces between each bodies in a matrix of size NxN 
std::vector<std::vector<double>> fixed_seq_forces(std::vector<Body>& bodies, size_t N){
    std::vector<std::vector<double>> results(N, std::vector<double>(N, 1));
    for (size_t i=0;i<N; i++){
        for (size_t j=0;j<N; j++){
            double f;
            if (i==j){
                f = 0.;
            }else if (results[i][j]!=0){
                double dx = bodies[i].x- bodies[j].x;
                double dy =  bodies[i].y- bodies[j].y;
                double mass_product = bodies[i].mass*bodies[j].mass;

                f = G * mass_product / (square(dx)+square(dy));
                results[i][j] = f;
                results[j][i] = f;
            }
        }
    }
    return results; 
}
