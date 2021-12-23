//
// Created by Guochenlin on 12/23/2021.
//

#include "main.h"
#include <iostream>
#include <stdio.h>
#include "FDM/Mesh.h"
#include <map>

int main(){
    Mesh mesh;
    std::map<std::string,std::string> Config =  mesh.ReadConfiguration("../Config.txt");
    mesh.GenerateMesh(Config);
    for(auto point:mesh.points){
        std::cout<<point<<std::endl;
    }





}