#ifndef GM_h
#define GM_h
#include "Dense"
#include <iostream>
#include <vector>

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

class GM
{
    public:
        const double PI  =3.141592653589793238463;
        int dataNum;
        double dt;
        double dur;
        std::string label;
        double PGA;
        double PGA_Scaled;
        double unit;
        double scaleFactor;
        bool isSpectrum;
        bool isScale;
        ArrayXd Ts;
        ArrayXd Sa;
        ArrayXd Sd;
        ArrayXd Sv;
        ArrayXd Sa_s;
        ArrayXd Sd_s;
        ArrayXd Sv_s;
        ArrayXd times;
        ArrayXd accel;
        ArrayXd accelScaled;
        GM(ArrayXd accel, double dt, double unit=9.81, std::string label="A single ground motion");
        ~GM();
        void Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
        vector<ArrayXXd> Response(double T, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
        void reScale(double scaleFactor);
        void scalePGA(double target);
        void saveTxt(std::string path);
private:
        vector<ArrayXXd> Newmark(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
        void update();
        void Scale();

};

#endif
