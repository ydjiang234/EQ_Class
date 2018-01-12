#include <iostream>
#include <vector>
#include "GM.h"
#include "FileLoadSaver.h"
#include "EigenArrayConvertor.h"
#include "GMSet.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

int main()
{
    FileLoadSaver *FLS = new FileLoadSaver();
    EigenArrayConvertor *Con = new EigenArrayConvertor();


    string path = "D:/Google Drive/Frames/GroundM/1-Chi-ChiTaiwan-04.dat";
    ArrayXXd accel = Con->OneDimArray(FLS->FileToDoubleArray(path));
    double dt = 0.001;
    double unit = 9.81;
    string label = "Test GM1";
    ArrayXd Ts = ArrayXd::LinSpaced(40, 0.01, 4.0);
    GM *gm1 = new GM(accel, dt, unit, label);

    /*
    gm1->Spectrum(Ts);
    string outpath = "C:/Users/yjiang/Desktop/tmp.txt";
    gm1->saveTxt(outpath);
    */
    GMSet *gmset = new GMSet();
    gmset->appendGM(*gm1);
    gmset->append(accel, 0.05, unit, "Test GM2");
    gmset->Spectrum(Ts);
    //std::cout<<gmset->GMs[0].label<<std::endl;
    //std::cout<<gmset->GMs[1].label<<std::endl;
    //std::cout<<gmset->GMNum<<std::endl;
    return 0;
}
