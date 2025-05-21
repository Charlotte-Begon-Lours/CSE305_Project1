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

};

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

int main(){
    std::vector<Body> bodies;
    
    //Body(double mass, double x, double y, double vX, double vY)
    Body earth = Body(5.97237e24,0, 0, 0, 0);
    earth.idBody = 0;
    Body moon = Body(7.342e22, 384400000.0, 0, 0, 0);
    moon.idBody = 1;
    bodies.push_back(earth);
    bodies.push_back(moon);
    size_t N = bodies.size();
    std::vector<std::vector<double>> results(N, std::vector<double>(N, 0));
    results = fixed_seq_forces(bodies, N);
    std::cout << "Force on object earth and moon : " << results[0][1] << " N" << std::endl;

    for (int i = 0; i < N; ++i) {
        std::cout << "Particle " << bodies[i].idBody  << ":\n";
        std::cout << "  Mass = " << bodies[i].mass << "\n";
        std::cout << "  Coordinates = (" << bodies[i].x << ", " << bodies[i].y << ")\n";
        std::cout << "--------------------------\n";
        std::cout << "Force on object earth and moon " << i << ": " << results[0][1] << " N" << std::endl;
    }

    return 0;
}