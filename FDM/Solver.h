//
// Created by Guochenlin on 12/26/2021.
//

#ifndef GSOLVER_SOLVER_H
#define GSOLVER_SOLVER_H

#include "Mesh.h"
#include "Init.h"
#include "Field.h"
#include <vector>
class Solver {

public:
    class ico{

    public:
        double dt;

        class projection: public ico{
        public:
            bool implicit=false;

            projection(bool implicit){
                this->implicit = implicit;
            }

            void predictor(Init &init,Mesh &mesh){
                Field::field Vxstar,Vystar;
                Vxstar.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
                Vystar.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
                Field::field & Ux = init.FieldList.find("Ux")->second;
                Field::field & Uy = init.FieldList.find("Uy")->second;

                for(int j = 0 ; j<mesh.Y_seed;j++){
                    for(int i=0; i<mesh.X_seed){
                        Vxstar.field[j][i] = Ux.field[j][i] +  * (Conv())
                    }
                }

            }


        };


        void Conv(Init &init,Mesh &mesh,std::string var,std::string scheme = "first_order_up_wind"){

            Field::field convx_field,convy_field;
            convx_field.name = var+"convx";
            convy_field.name = var+"convy";

            convx_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
            convy_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
            if(scheme=="first_order_up_wind"){

                Field::field work_field = init.FieldList.find(var)->second;
                Field::field Ux = init.FieldList.find("Ux")->second;
                Field::field Uy = init.FieldList.find("Uy")->second;

                for(int j=1;j<mesh.Y_seed-1;j++) {
                    for (int i = 1; i < mesh.X_seed - 1; i++) {
                        /*
                         *       idn
                         *   idw|id|ide
                         *       ids
                         *
                         */
                        int id = mesh.IJKtoID(i,j);
                        int ide = mesh.IJKtoID(i+1,j);
                        int idw = mesh.IJKtoID(i-1,j);
                        int idn = mesh.IJKtoID(i,j+1);
                        int ids = mesh.IJKtoID(i,j-1);
                        convx_field.field[j][i] = Ux.field[j][i] * (work_field.field[j][i+1] - work_field.field[j][i] )/(mesh.points[ide].location[0] - mesh.points[id].location[0] );
                        convy_field.field[j][i] = Uy.field[j][i] * (work_field.field[j+1][i] - work_field.field[j-1][i] )/(mesh.points[idn].location[0] - mesh.points[ids].location[0] );
                    }
                }
                init.FieldList.insert(std::pair<std::string,Field::field>(var+"convx",convx_field));
                init.FieldList.insert(std::pair<std::string,Field::field>(var+"convy",convy_field));

            }else{
                std::cout<<"The scheme is not implemented yet!";
            }

        }

        void Diff(Init &init,Mesh &mesh,std::string var,std::string scheme = "second_order_central_diff"){

            Field::field diffx_field,diffy_field;
            diffx_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
            diffy_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
            if(scheme=="second_order_central_diff"){

                Field::field work_field = init.FieldList.find(var)->second;
                for(int j=1;j<mesh.Y_seed-1;j++) {
                    for (int i = 1; i < mesh.X_seed - 1; i++) {
                        /*
                         *       idn
                         *   idw|id|ide
                         *       ids
                         *
                         */
                        int id = mesh.IJKtoID(i,j);
                        int ide = mesh.IJKtoID(i+1,j);
                        int idw = mesh.IJKtoID(i-1,j);
                        int idn = mesh.IJKtoID(i,j+1);
                        int ids = mesh.IJKtoID(i,j-1);
                        diffx_field.field[j][i] = (work_field.field[j][i-1] - 2 * work_field.field[j][i] + work_field.field[j][i+1] ) / pow(mesh.points[ide].location[0] - mesh.points[idw].location[0],2);
                        diffy_field.field[j][i] = (work_field.field[j+1][i] - 2 * work_field.field[j][i] + work_field.field[j-1][i] ) / pow(mesh.points[idn].location[0] - mesh.points[ids].location[0],2);
                    }
                }
                init.FieldList.insert(std::pair<std::string,Field::field>(var+"diffx",diffx_field));
                init.FieldList.insert(std::pair<std::string,Field::field>(var+"diffy",diffy_field));

            }else{
                std::cout<<"The scheme is not implemented yet!";
            }

        }


        //time loops
        void timeloops(Mesh &mesh) {

            int totalTimesteps = std::stoi(mesh.FindConfig(mesh.keyValueMap,"totalTimesteps"));
            for (int timesteps = 0; timesteps < totalTimesteps; timesteps++) {


            }


        }


    };


};


#endif //GSOLVER_SOLVER_H
