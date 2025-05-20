#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double G = 6.67430e-11;
 
class Body{
    int idBody;
    double mass;
    double x, y;
    double vX, vY;
    double aX, aY;
    Body(double mass, double x, double y, double vX, double vY)
    : mass(mass), x(x), y(y), vX(vX), vY(vY), aX(0), aY(0){}
};

int main(){
    return 0;
}