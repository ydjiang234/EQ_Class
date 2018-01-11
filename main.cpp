#include <iostream>
#include <vector>
#include "GM.h"
#include "FileLoader.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using std::vector;

int main()
{
    /*
    ArrayXd accel;
    accel = ArrayXd(3,1);
    accel<<1,2,3;
    double dt = 0.05;
    GM *gm1 = new GM(accel, dt, 9.81, "90");
    std::cout<<gm1->accel<<std::endl;
    delete gm1;
    */
    char path[100] = "D:/Google Drive/Python_Scripts/Test/ABAQUS/1_FE.txt";
    FileLoader *FL = new FileLoader();

    std::cout<<FL->FileToLongArray(path)[50][1]-100<<std::endl;
    return 0;
}
