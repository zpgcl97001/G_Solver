//
// Created by Guochenlin on 12/26/2021.
//

#ifndef GSOLVER_POST_H
#define GSOLVER_POST_H
#include "Field.h"
#include "Mesh.h"
#include "Init.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>

class Post {

public:

    void ExportTecplot(Init &init,Mesh &mesh,std::string filename="result.dat"){

        if(filename != "0"){
            filename += ".dat";
        }

        std::ofstream outfile;
        outfile.open(filename,std::ios::out);
        outfile<<std::setprecision(14);
        outfile.setf(std::ios::scientific);



        outfile << "VARIABLES =\"x\" \n"
                << "\"y\" \n";

        for(auto field_pair : init.FieldList){
            outfile<<"\" " + field_pair.first + "\"\n";
        }
        outfile<<"ZONE T = \"Rank" << 1<<"\" \n"
               <<"I= "<<std::to_string(mesh.X_seed)<<"\n"
               <<"J= "<<std::to_string(mesh.Y_seed)<<"\n";
        outfile<<"DATAPACKING = POINT\n";

        for(int i = 0 ; i<mesh.Y_seed;i++){
            for(int j = 0;j<mesh.X_seed;j++){
                outfile<<mesh.points[mesh.IJKtoID(j,i)].location[0]<<" "<<mesh.points[mesh.IJKtoID(j,i)].location[1]<<" ";
                std::cout<<j<<","<<i<<","<<mesh.IJKtoID(j,i)<<std::endl;
                for(auto field_pair : init.FieldList){
                    outfile<<field_pair.second.field[i][j]<<" ";
                }
                outfile<<"\n";
            }
        }
        outfile.unsetf(std::ios::scientific);
        outfile << std::setprecision(6);
        outfile.close();

        }



};


#endif //GSOLVER_POST_H
