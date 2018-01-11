#include "GM.h"
#include "Dense"
#include <iostream>

using Eigen::ArrayXd;
using Eigen::ArrayXXd;

GM::GM(ArrayXd accel, double dt, double unit, std::string label)
{
    this->accel = accel;
    this->dt = dt;
    this->unit = unit;
    this->label = label;
}
GM::~GM()
{
    std::cout<<this->label;
    std::cout<<" is deleted"<<std::endl;
}

void GM::Spectrum(ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    int a = 1;
}

ArrayXXd GM::Newmark(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    ArrayXd a(5,1);
    return a;
}

ArrayXXd GM::response(double T, double dampingRatio, double beta, double gamma)
{
    ArrayXd a(5,1);
    return a;
}
