//
// Created by Guochenlin on 12/26/2021.
//

#ifndef GSOLVER_SOLVER_H
#define GSOLVER_SOLVER_H

#include "Mesh.h"
#include "Init.h"
#include "Field.h"
#include <vector>
#include <armadillo>
class Solver {

public:


};

class ico: public Solver{

public:
    double dt,mu,rho;


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
                    convx_field.field[j][i] = Ux.field[j][i] * (work_field.field[j][i] - work_field.field[j][i-1] )/(mesh.points[ide].location[0] - mesh.points[id].location[0] );
                    convy_field.field[j][i] = Uy.field[j][i] * (work_field.field[j+1][i] - work_field.field[j-1][i] )/(mesh.points[idn].location[1] - mesh.points[ids].location[1] );
                }
            }
            init.FieldList.insert(std::pair<std::string,Field::field>(var+"convUx",convx_field));
            init.FieldList.insert(std::pair<std::string,Field::field>(var+"convUy",convy_field));

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
                    diffy_field.field[j][i] = (work_field.field[j+1][i] - 2 * work_field.field[j][i] + work_field.field[j-1][i] ) / pow(mesh.points[idn].location[1] - mesh.points[ids].location[1],2);
                }
            }
            init.FieldList.insert(std::pair<std::string,Field::field>(var+"diffUx",diffx_field));
            init.FieldList.insert(std::pair<std::string,Field::field>(var+"diffUy",diffy_field));

        }else{
            std::cout<<"The scheme is not implemented yet!";
        }

    }

    void Div(Init &init,Mesh &mesh,std::string Ux,std::string Uy,std::string scheme = "second_order_central"){

        Field::field div_field;
        div_field.name = "div"+Ux[0];
        div_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));

        Field::field work_fieldx = init.FieldList.find(Ux)->second;
        Field::field work_fieldy = init.FieldList.find(Uy)->second;

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
                div_field.field[j][i] =  (work_fieldx.field[j][i+1]-work_fieldx.field[j][i-1])/(mesh.points[ide].location[0]-mesh.points[idw].location[0]);
                div_field.field[j][i] += (work_fieldy.field[j+1][i]-work_fieldy.field[j-1][i])/(mesh.points[idn].location[1]-mesh.points[ids].location[1]);
            }
        }



    }


    //time loops
    void timeloops(Mesh &mesh) {

        int totalTimesteps = std::stoi(mesh.FindConfig(mesh.keyValueMap,"totalTimesteps"));
        for (int timesteps = 0; timesteps < totalTimesteps; timesteps++) {

        }
    }

};

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
        std::vector<std::vector<double>> & Ux = init.FieldList.find("Ux")->second.field;
        std::vector<std::vector<double>> & Uy = init.FieldList.find("Uy")->second.field;
        this->Conv(init,mesh,"Ux");
        this->Conv(init,mesh,"Uy");
        this->Diff(init,mesh,"Ux");
        this->Diff(init,mesh,"Uy");
        std::vector<std::vector<double>> &convUx = init.FieldList.find("convUx")->second.field;
        std::vector<std::vector<double>> &convUy = init.FieldList.find("convUy")->second.field;
        std::vector<std::vector<double>> &diffUx = init.FieldList.find("convUy")->second.field;
        std::vector<std::vector<double>> &diffUy = init.FieldList.find("convUy")->second.field;



        for(int j = 0 ; j<mesh.Y_seed;j++){
            for(int i=0; i<mesh.X_seed;i++){
                Vxstar.field[j][i] = Ux[j][i] + dt * (-convUx[j][i] + mu * diffUx[j][i]);
                Vystar.field[j][i] = Uy[j][i] + dt * (-convUx[j][i] + mu * diffUx[j][i]);
            }
        }
    }

    void Possion(Init &init,Mesh &mesh){

        Div(init,mesh,"Vxstar","Vystar");
        std::vector<std::vector<double>> DivUstar = init.FieldList.find("divU")->second.field;
        std::vector<std::vector<double>> RHS(mesh.Y_seed,std::vector<double> (mesh.X_seed));

        for(int j=0;j<mesh.Y_seed;j++){
            for(int i=0;i<mesh.X_seed;i++)
            {
                RHS[j][i] = DivUstar[j][i] * rho/dt;
            }
        }

        std::vector<std::vector<double>> LaplaceP(mesh.Y_seed*mesh.X_seed,std::vector<double> (mesh.Y_seed*mesh.X_seed,0));

        double X_length = std::stof(mesh.FindConfig(mesh.keyValueMap,"X_length"));
        double Y_length = std::stof(mesh.FindConfig(mesh.keyValueMap,"Y_length"));
        double X_seed = std::stof(mesh.FindConfig(mesh.keyValueMap,"X_seed"));
        double Y_seed = std::stof(mesh.FindConfig(mesh.keyValueMap,"Y_seed"));
        double dx = X_length/X_seed;
        double dy = Y_length/Y_seed;

        for(int id = 0; id<mesh.X_seed*mesh.Y_seed;id++){

            std::vector<int> IJK = mesh.IDtoIJK(id);
            if(IJK[0]==0){

            }
            else if(IJK[0] == mesh.X_seed - 1){

            }else if(IJK[1] == 0){

            }else if(IJK[1] == mesh.Y_seed - 1){

            }else{
                /*
                 *                    Pn
                 *         Pw    |    P    |   Pe
                 *                    Ps
                 */

                int ide = mesh.IJKtoID(IJK[0]+1,IJK[1],0);
                int idw = mesh.IJKtoID(IJK[0]-1,IJK[1],0);
                int idn = mesh.IJKtoID(IJK[0],IJK[1]+1,0);
                int ids = mesh.IJKtoID(IJK[0],IJK[1]-1,0);


                LaplaceP[id][id] += -2/pow(dx,2)  -2/pow(dy,2);
                LaplaceP[id][idw] += 1/pow(dx,2);
                LaplaceP[id][ide] += 1/pow(dy,2);
                LaplaceP[id][idn] += 1/pow(dy,2);
                LaplaceP[id][ids] += 1/pow(dy,2);

            }

        }

        arma::mat A = arma::conv_to<arma::mat>::from(LaplaceP);
        arma::mat b = arma::conv_to<arma::mat>::from(RHS);
        b = b.t();
        arma::mat x = arma::solve(A,b);



    }


};

#endif //GSOLVER_SOLVER_H
