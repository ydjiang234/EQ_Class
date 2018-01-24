#include "GM.h"
#include "Dense"
#include <iostream>
#include <cmath>
#include <vector>
#include "EArrayIO.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using std::pow;
using std::vector;

GM::GM(ArrayXd accel, double dt, double unit, std::string label)
{
    this->accel = accel;
    this->dt = dt;
    this->unit = unit;
    this->label = label;
    this->dataNum = accel.rows();
    this->scaleFactor = 1.0;
    this->isSpectrum = false;
    this->isScale = false;
    this->times = this->dt * ArrayXd::LinSpaced(this->dataNum, 0, this->dataNum);
    this->Scale();
}
GM::~GM()
{
    //std::cout<<this->label;
    //std::cout<<" is deleted"<<std::endl;
}

void GM::reScale(double scaleFactor)
{
    this->scaleFactor = scaleFactor;
    this->Scale();
}

void GM::Scale()
{
    this->accelScaled = this->accel * this->scaleFactor;
    this->PGA = this->accel.maxCoeff();
    this->PGA_Scaled = this->accelScaled.maxCoeff();
    if (this->isSpectrum)
    {
        this->Sa_s = this->Sa * this->scaleFactor;
        this->Sv_s = this->Sv * this->scaleFactor;
        this->Sd_s = this->Sd * this->scaleFactor;
    }

    this->isScale = true;
}

void GM::SpectrumFull(ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    this->Ts = Ts;
    vector<ArrayXXd> input = this->NewmarkFull(this->dt, this->accel*this->unit, this->Ts, dampingRatio, beta, gamma);
    input[2] = input[2] / this->unit;
    input[3] = input[3] / this->unit;
    this->Sa = input[3].abs().colwise().maxCoeff();
    this->Sv = input[0].abs().colwise().maxCoeff();
    this->Sd = input[1].abs().colwise().maxCoeff();
    this->isSpectrum = true;
    this->Scale();
}

void GM::Spectrum(ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    this->Ts = Ts;
    vector<ArrayXd> input = this->Newmark(this->dt, this->accel*this->unit, this->Ts, dampingRatio, beta, gamma);

    this->Sa = input[0] / this->unit;
    this->Sv = input[1];
    this->Sd = input[2];
    this->isSpectrum = true;
    this->Scale();

}

void GM::scalePGA(double target)
{
    this->reScale(target/this->PGA);
}

vector<ArrayXXd> GM::NewmarkFull(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    int dataNum = accelG.rows();
    int Tnum = Ts.rows();
    double m = 1.0;//mass
    ArrayXd k = 4.0 * pow(this->PI,2) * m / Ts.pow(2);
    ArrayXd cc = 2.0 * (k * m).sqrt();
    ArrayXd c = cc * dampingRatio;
    ArrayXd meq = m + c * gamma * dt + k * beta * pow(dt,2);
    ArrayXXd a_s = ArrayXXd::Zero(dataNum, Tnum);
    ArrayXXd v_s = ArrayXXd::Zero(dataNum, Tnum);
    ArrayXXd d_s = ArrayXXd::Zero(dataNum, Tnum);
    for (int i=1; i<dataNum; i++)
    {
        double dag = accelG(i) - accelG(i-1);
        ArrayXd a1 = a_s.row(i-1);
        ArrayXd v1 = v_s.row(i-1);
        ArrayXd d1 = d_s.row(i-1);
        ArrayXd dd = v1 * dt + a1 * pow(dt,2) / 2.0;
        ArrayXd dv = a1 * dt;
        ArrayXd da = (-1.0 * m * dag - c * dv - k * dd) / meq;
        a_s.row(i) = a1 + da;
        v_s.row(i) = v1 + dv + da * gamma * dt;
        d_s.row(i) = d1 + dd + da * beta * pow(dt,2);
    }
    ArrayXXd a_s_abs = a_s + (accelG.matrix() * MatrixXd::Ones(1,Tnum)).array();
    vector<ArrayXXd> out = {v_s, d_s, a_s, a_s_abs};
    return out;
}


vector<ArrayXd> GM::Newmark(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    int dataNum = accelG.rows();
    int Tnum = Ts.rows();
    double m = 1.0;//mass
    ArrayXd k = 4.0 * pow(this->PI,2) * m / Ts.pow(2);
    ArrayXd cc = 2.0 * (k * m).sqrt();
    ArrayXd c = cc * dampingRatio;
    ArrayXd meq = m + c * gamma * dt + k * beta * pow(dt,2);
    ArrayXd a1(Tnum), v1(Tnum), d1(Tnum);
    ArrayXd a2(Tnum), v2(Tnum), d2(Tnum);
    ArrayXd a_abs1(Tnum), a_abs2(Tnum);
    ArrayXd Sa(Tnum), Sv(Tnum), Sd(Tnum);
    a1 = ArrayXd::Zero(Tnum);
    v1 = ArrayXd::Zero(Tnum);
    d1 = ArrayXd::Zero(Tnum);
    a_abs1 = a1 + accelG(0);
    for (int i=1; i<dataNum; i++)
    {
        double dag = accelG(i) - accelG(i-1);
        ArrayXd dd = v1 * dt + a1 * pow(dt,2) / 2.0;
        ArrayXd dv = a1 * dt;
        ArrayXd da = (-1.0 * m * dag - c * dv - k * dd) / meq;
        a2 = a1 + da;
        v2 = v1 + dv + da * gamma * dt;
        d2 = d1 + dd + da * beta * pow(dt,2);
        a_abs2 = a2 + accelG(i);
        if (i==1)
        {
            Sa = this->absmax(a_abs1,a_abs2);
            Sv = this->absmax(v1,v2);
            Sd = this->absmax(d1,d2);
        } else {
            Sa = this->absmax(Sa,a_abs2);
            Sv = this->absmax(Sv,v2);
            Sd = this->absmax(Sd,d2);
        }
        a1 = a2;
        v1 = v2;
        d1 = d2;
    }
    vector<ArrayXd> out = {Sa, Sv, Sd};
    return out;
}

ArrayXd GM::absmax(ArrayXd a, ArrayXd b)
{
    return a.abs().cwiseMax(b.abs());
}


vector<ArrayXXd> GM::Response(double T, double dampingRatio, double beta, double gamma)
{
    ArrayXd curT(1);
    curT(0) = T;
    vector<ArrayXXd> out = this->NewmarkFull(this->dt, this->accel*this->unit, curT, dampingRatio, beta, gamma);
    out[2] = out[2] / this->unit;
    out[3] = out[3] / this->unit;
    return out;
}

void GM::savetxtRespone(std::string path)
{
    if (this->isSpectrum)
    {
        EArrayIO *EAIO = new EArrayIO();

        ArrayXXd out(this->Ts.rows(), 4);
        out.col(0) = this->Ts;
        out.col(1) = this->Sa_s;
        out.col(2) = this->Sv_s;
        out.col(3) = this->Sd_s;
        EAIO->savetxt(out, path);
    }
}

void GM::savetxtAccel(std::string path)
{
    EArrayIO *EAIO = new EArrayIO();
    EAIO->savetxt(this->accelScaled, path);
}


//backup
/*
vector<ArrayXXd> GM::Newmark(double dt, ArrayXd accelG, ArrayXd Ts, double dampingRatio, double beta, double gamma)
{
    int dataNum = accelG.rows();
    int Tnum = Ts.rows();
    double m = 1.0;//mass
    ArrayXd k = 4.0 * pow(this->PI,2) * m / Ts.pow(2);
    ArrayXd cc = 2.0 * (k * m).sqrt();
    ArrayXd c = cc * dampingRatio;
    ArrayXd meq = m + c * gamma * dt + k * beta * pow(dt,2);
    ArrayXXd a_s = ArrayXXd::Zero(dataNum, Tnum);
    ArrayXXd v_s = ArrayXXd::Zero(dataNum, Tnum);
    ArrayXXd d_s = ArrayXXd::Zero(dataNum, Tnum);
    for (int i=1; i<dataNum; i++)
    {
        double dag = accelG(i) - accelG(i-1);
        ArrayXd a1 = a_s.row(i-1);
        ArrayXd v1 = v_s.row(i-1);
        ArrayXd d1 = d_s.row(i-1);
        ArrayXd dd = v1 * dt + a1 * pow(dt,2) / 2.0;
        ArrayXd dv = a1 * dt;
        ArrayXd da = (-1.0 * m * dag - c * dv - k * dd) / meq;
        a_s.row(i) = a1 + da;
        v_s.row(i) = v1 + dv + da * gamma * dt;
        d_s.row(i) = d1 + dd + da * beta * pow(dt,2);
    }
    ArrayXXd a_s_abs = a_s + (accelG.matrix() * MatrixXd::Ones(1,Tnum)).array();
    vector<ArrayXXd> out = {v_s, d_s, a_s, a_s_abs};
    return out;
}
*/
