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

    double unit = 9.81;
    ArrayXd Ts = ArrayXd::LinSpaced(400, 0.01, 4.0);
    
    std::string GM_Path = "D:/Google Drive/Frames/GroundM/Original";
    std::string GM_ID_Path = GM_Path + "/GM_ID.txt";
    std::string GM_dt_Path = GM_Path + "/GM_dt.txt";
    std::string GM_scaleFactor_Path = GM_Path + "/Scale_factor.txt";

    //Load the GM_id
    vector<vector<std::string>> GM_ID = FLS->FileToStringArray(GM_ID_Path);
    //Load GM_dt
    vector<vector<double>> GM_dt = FLS->FileToDoubleArray(GM_dt_Path);
    //Load scale factor
    vector<vector<double>> temp = FLS->FileToDoubleArray(GM_scaleFactor_Path);
    vector<double> scaleFactors;
    for (int i=0;i<temp.size();i++)
    {
        scaleFactors.push_back(temp[i][0]);
    }
       
    GMSet *gmset = new GMSet();
    for (int i=0; i<GM_ID.size();i++)
    {
        ArrayXd accel = Con->OneDimArray(FLS->FileToDoubleArray(GM_Path + '/' + GM_ID[i][0]));
        gmset->append(accel, GM_dt[i][0], unit, GM_ID[i][0]);
        std::cout<<GM_ID[i][0]<<" is Loaded."<<std::endl;
    }


    //gmset->Spectrum(Ts);

    //gmset->savetxt(GM_Path + "/Set");
    //std::cout<<gmset->GMs[0].label<<std::endl;
    //std::cout<<gmset->GMs[1].label<<std::endl;
    //std::cout<<gmset->Sa<<std::endl;
    //gmset->loadtxt(GM_Path + "/Set");
    //gmset->reScale(scaleFactors);
    //gmset->savetxt(GM_Path + "/Set_Scaled");

    return 0;
}
