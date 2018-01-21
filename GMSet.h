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
    bool isSpectrum;
    vector<double> dts;
    ArrayXd Ts;
    ArrayXXd Sa;
    ArrayXXd Sv;
    ArrayXXd Sd;
    ArrayXd SaMean;
    ArrayXd SdMean;
    ArrayXd SvMean;
    vector<std::string> labels;
    vector<double> scaleFactors;
    void appendGM(GM gm);
    void append(ArrayXd accel, double dt, double unit=9.81, std::string label="A single GroundMotion");
    void Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
    void update();
    void savetxt(std::string path);
    void loadtxt(std::string path);
    void scale();
    void reScale(vector<double> input);

    GMSet();
    ~GMSet();
};

#endif // GMSET_H
