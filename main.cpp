//
// Created by Guochenlin on 12/23/2021.
//

#include "main.h"
#include <iostream>
#include <stdio.h>
#include "FDM/Mesh.h"
#include "FDM/Field.h"
#include "FDM/Init.h"
#include "FDM/Post.h"
#include <map>
#include <chrono>

int main() {

    auto start_time = std::chrono::steady_clock::now();
    Mesh mesh;
    Field fields;
    Init init;
    Post post;
    std::map<std::string, std::string> Config = mesh.ReadConfiguration("../Config.txt");
    mesh.GenerateMesh(Config);
    init.MammalInit(mesh);
    post.ExportTecplot(init,mesh);

    auto end_time =  std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    std::cout<<"Init Finished! Duration = "<<duration.count() <<" s"<<std::endl;
    return 0;

}
