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
    ArrayXd Ts = ArrayXd::LinSpaced(20, 0.01, 4.0);
    
    std::string GM_ID_Path = "./Original/GM_ID.txt";
    std::string GM_dt_Path = "./Original/GM_dt.txt";
    std::string GM_scaleFactor_Path = "./Original/Scale_factor.txt";

    //Load the GM_id
    vector<vector<std::string>> GM_ID = FLS->FileToStringArray(GM_ID_Path);
    //Load GM_dt
    vector<vector<double>> GM_dt = FLS->FileToDoubleArray(GM_dt_Path);
    //Load scale factor
    vector<vector<double>> scaleFactors = FLS->FileToDoubleArray(GM_scaleFactor_Path);
    
       
    GMSet *gmset = new GMSet();
    for (int i=0; i<GM_ID.size();i++)
    {
        ArrayXd accel = Con->OneDimArray(FLS->FileToDoubleArray("./Original/" + GM_ID[i][0]));
        gmset->append(accel, GM_dt[i][0], unit, GM_ID[i][0]);
    }


    //gmset->Spectrum(Ts);
    //gmset->savetxt("Set");
    //std::cout<<gmset->GMs[0].label<<std::endl;
    //std::cout<<gmset->GMs[1].label<<std::endl;
    //std::cout<<gmset->Sa<<std::endl;
    gmset->loadtxt("Set");
    std::cout<<gmset->Sa<<std::endl;
    return 0;
}
