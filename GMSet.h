#ifndef GMSET_H
#define GMSET_H
#include <iostream>
#include <vector>
#include "GM.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

class GMSet
{
public:
    vector<GM> GMs;
    int GMNum;
    vector<double> dts;
    ArrayXd Ts;
    vector<std::string> labels;
    vector<double> scaleFactors;
    void appendGM(GM gm);
    void append(ArrayXd accel, double dt, double unit=9.81, std::string label="A single GroundMotion");
    void Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
    void savetxt(std::string path);
    void loadtxt(std::string path);

    GMSet();
    ~GMSet();
};

#endif // GMSET_H
