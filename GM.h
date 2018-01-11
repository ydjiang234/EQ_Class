#ifndef GM_h
#define GM_h
#include "Dense"
#include <iostream>

using Eigen::ArrayXd;
using Eigen::ArrayXXd;

class GM
{
    public:
        int dataNum;
        double dt;
        double dur;
        std::string label;
        double PGA;
        double unit;
        double scaleFactor;
        bool ifSpectrum;
        ArrayXd times;
        ArrayXd accel;
        ArrayXd accleScaled;
        GM(ArrayXd accel, double dt, double unit=9.81, std::string label="A single ground motion");
        ~GM();
        void Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);

    private:
        ArrayXXd Newmark(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
        ArrayXXd response(double T, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
};

#endif
