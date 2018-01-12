#include "GMSet.h"
#include <iostream>
#include <vector>
#include "GM.h"
#include "FileLoadSaver.h"
#include "EigenArrayConvertor.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

GMSet::GMSet()
{
    this->GMNum = 0;
}

GMSet::~GMSet() {}


void GMSet::appendGM(GM gm)
{
    this->GMs.push_back(gm);
    this->dts.push_back(gm.dt);
    this->labels.push_back(gm.label);
    this->scaleFactors.push_back(gm.scaleFactor);
    this->GMNum += 1;
}

void GMSet::append(ArrayXd accel, double dt, double unit, std::string label)
{
    GM curGM = GM(accel, dt, unit, label);
    this->appendGM(curGM);
}

void GMSet::Spectrum(ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    for (int i=0; i<this->GMNum; i++)
    {
        this->GMs[i].Spectrum(Ts, dampingRatio, beta, gamma);
    }
}

/*
void GMSet::Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
void GMSet::savetxt(std::string path);
void GMSet::loadtxt(std::string path);
*/

