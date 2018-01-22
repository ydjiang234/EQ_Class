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
    this->isSpectrum = false;
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
    this->Ts = Ts;
    for (int i=0; i<this->GMNum; i++)
    {
        this->GMs[i].Spectrum(Ts, dampingRatio, beta, gamma);
    }
    this->isSpectrum = true;
    this->update();
}

void GMSet::update()
{
    this->Sa = ArrayXXd::Zero(this->Ts.rows(), this->GMNum);
    this->Sv = ArrayXXd::Zero(this->Ts.rows(), this->GMNum);
    this->Sd = ArrayXXd::Zero(this->Ts.rows(), this->GMNum);
    this->SaMean= ArrayXd::Zero(this->Ts.rows());
    this->SvMean= ArrayXd::Zero(this->Ts.rows());
    this->SdMean= ArrayXd::Zero(this->Ts.rows());
    for (int i=0; i<this->GMNum; i++)
    {
        this->Sa.col(i) = this->GMs[i].Sa_s;
        this->Sv.col(i) = this->GMs[i].Sv_s;
        this->Sd.col(i) = this->GMs[i].Sd_s;
    }
        this->SaMean = this->Sa.rowwise().mean();
        this->SvMean = this->Sv.rowwise().mean();
        this->SdMean = this->Sd.rowwise().mean();
}

void GMSet::savetxt(std::string path)
{
    ArrayXXd outSa = ArrayXXd::Zero(this->Ts.rows(), this->GMNum+2);
    ArrayXXd outSv = ArrayXXd::Zero(this->Ts.rows(), this->GMNum+2);
    ArrayXXd outSd = ArrayXXd::Zero(this->Ts.rows(), this->GMNum+2);
    outSa << this->Ts, this->SaMean, this->Sa;
    outSv << this->Ts, this->SvMean, this->Sv;
    outSd << this->Ts, this->SdMean, this->Sd;
    FileLoadSaver *FLS = new FileLoadSaver();
    EigenArrayConvertor *Con = new EigenArrayConvertor();
    FLS->DoubleArrayToFile(Con->FromTwoDimArray(outSa), path + "_Sa.txt", ' ');
    FLS->DoubleArrayToFile(Con->FromTwoDimArray(outSv), path + "_Sv.txt", ' ');
    FLS->DoubleArrayToFile(Con->FromTwoDimArray(outSd), path + "_Sd.txt", ' ');
    std::cout<<"OK"<<std::endl;
}

void GMSet::loadtxt(std::string path)
{
    FileLoadSaver *FLS = new FileLoadSaver();
    EigenArrayConvertor *Con = new EigenArrayConvertor();
    ArrayXXd inSa = Con->TwoDimArray(FLS->FileToDoubleArray(path + "_Sa.txt"));
    ArrayXXd inSv = Con->TwoDimArray(FLS->FileToDoubleArray(path + "_Sv.txt"));
    ArrayXXd inSd = Con->TwoDimArray(FLS->FileToDoubleArray(path + "_Sd.txt"));
    this->Ts = inSa.col(0);
    for (int i=0; i<this->GMNum; i++)
    {
        this->GMs[i].Sa = inSa.col(i+2);
        this->GMs[i].Sv = inSv.col(i+2);
        this->GMs[i].Sd = inSd.col(i+2);
        this->GMs[i].isSpectrum = true;
        this->GMs[i].Scale();
    }
    this->isSpectrum = true;
    this->update();
}

void GMSet::scale()
{
    for (int i=0; i<this->GMNum; i++)
    {
        this->GMs[i].reScale(this->scaleFactors[i]);
    }
    this->update();

}
void GMSet::reScale(vector<double> input)
{
    this->scaleFactors = input;
    this->scale();
}

/*
void GMSet::Spectrum(ArrayXd Ts, double dampingRatio=0.05, double beta=0.25, double gamma=0.5);
void GMSet::savetxt(std::string path);
void GMSet::loadtxt(std::string path);
*/

