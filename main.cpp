#include <iostream>
#include <vector>
#include "EArrayIO.h"
#include "VectorIO.h"
#include "GM.h"
#include "GMSet.h"
#include <chrono>


using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

int main()
{

    EArrayIO *EAIO = new EArrayIO();
    VectorIO *VIO = new VectorIO();

    double unit = 9.81;
    ArrayXd Ts = ArrayXd::LinSpaced(400, 0.01, 4.0);
    

    //Test GM
    /*
    std::string path = "D:/Google Drive/Frames/GroundM/Original/1-Chi-ChiTaiwan-04.dat";
    ArrayXd accel = EAIO->loadtxt(path);
    GM *gm1 = new GM(accel, 0.001, unit);
    //gm1->reScale(50.0);
    auto start = std::chrono::high_resolution_clock::now();
    gm1->Spectrum(Ts);
    auto finish = std::chrono::high_resolution_clock::now();
    //gm1->savetxtRespone("C:/Users/yjiang/Desktop/saved.txt");

    std::chrono::duration<double> elapsed = finish - start;
    std::cout<<elapsed.count()<<std::endl;
    //std::cout<<gm1->Sv<<std::endl;
    */


    //Test GMSet
    std::string GM_Path = "D:/Google Drive/Frames/GroundM/Original";
    std::string GM_ID_Path = GM_Path + "/GM_ID.txt";
    std::string GM_dt_Path = GM_Path + "/GM_dt.txt";
    std::string GM_scaleFactor_Path = GM_Path + "/Scale_factor.txt";

    //Load the GM_id
    vector<std::string> GM_ID = VIO->loadline(GM_ID_Path);
    //Load GM_dt
    ArrayXd GM_dt = EAIO->loadtxt(GM_dt_Path);
    //Load scale factor
    ArrayXd scaleFactors = EAIO->loadtxt(GM_scaleFactor_Path);


       
    GMSet *gmset = new GMSet();

    for (int i=0; i<GM_ID.size();i++)
    {
        ArrayXd accel = EAIO->loadtxt(GM_Path + "/" + GM_ID[i]);
        gmset->append(accel, GM_dt(i), unit, GM_ID[i]);
        //std::cout<<GM_ID[i]<<" is Loaded."<<std::endl;
    }


    gmset->Spectrum(Ts);

    gmset->savetxt(GM_Path + "/Set");
    //std::cout<<gmset->GMs[0].label<<std::endl;
    //std::cout<<gmset->GMs[1].label<<std::endl;
    //std::cout<<gmset->Sa<<std::endl;
    //gmset->loadtxt(GM_Path + "/Set");
    gmset->reScale(scaleFactors);
    gmset->savetxt(GM_Path + "/Set_Scaled");

    return 0;
}
