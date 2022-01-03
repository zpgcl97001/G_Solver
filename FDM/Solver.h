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
    double dt=1e-3,mu=0.01,rho=1;


    void Conv(Init &init,Mesh &mesh,std::string var,std::string scheme = "first_order_up_wind"){

        Field::field convx_field,convy_field;
        convx_field.name = var+"convx";
        convy_field.name = var+"convy";

        convx_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        convy_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        if(scheme=="first_order_up_wind"){

            Field::field &work_field = init.FieldList.find(var)->second;
            Field::field &Ux = init.FieldList.find("Ux")->second;
            Field::field &Uy = init.FieldList.find("Uy")->second;

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
                    convx_field.field[j][i] = Ux.field[j][i] * (work_field.field[j][i] - work_field.field[j][i-1] )/(mesh.points[id].location[0] - mesh.points[idw].location[0] );
                    convy_field.field[j][i] = Uy.field[j][i] * (work_field.field[j+1][i] - work_field.field[j-1][i] )/(mesh.points[idn].location[1] - mesh.points[ids].location[1] );
                }
            }
            init.FieldList.insert(std::pair<std::string,Field::field>("conv"+var,convx_field));
            init.FieldList.insert(std::pair<std::string,Field::field>("conv"+var,convy_field));

        }else{
            std::cout<<"The scheme is not implemented yet!";
        }

    }

    void Diff(Init &init,Mesh &mesh,std::string var,std::string scheme = "second_order_central_diff"){

        Field::field diffx_field,diffy_field;
        diffx_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        diffy_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        if(scheme=="second_order_central_diff"){

            Field::field &work_field = init.FieldList.find(var)->second;
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
            init.FieldList.insert(std::pair<std::string,Field::field>("diff"+var,diffx_field));
            init.FieldList.insert(std::pair<std::string,Field::field>("diff"+var,diffy_field));

        }else{
            std::cout<<"The scheme is not implemented yet!";
        }

    }

    void Div(Init &init,Mesh &mesh,std::string Ux,std::string Uy,std::string scheme = "second_order_central"){

        Field::field div_field;
        div_field.name = "divU";
        div_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));

        Field::field &work_fieldx = init.FieldList.find(Ux)->second;
        Field::field &work_fieldy = init.FieldList.find(Uy)->second;

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

        init.FieldList.insert(std::pair<std::string,Field::field>("divU",div_field));


    }


    void Grad(Init &init,Mesh &mesh,std::string Scalar,std::string scheme = "second_order_central"){

        Field::field gradx_field,grady_field;
        gradx_field.name = "gradx";
        grady_field.name = "grady";
        gradx_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        grady_field.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));

        Field::field &work_field = init.FieldList.find(Scalar)->second;

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
                gradx_field.field[j][i] = (work_field.field[j][i+1]-work_field.field[j][i-1])/(mesh.points[ide].location[0] - mesh.points[idw].location[0]);
                grady_field.field[j][i] = (work_field.field[j+1][i]-work_field.field[j-1][i])/(mesh.points[idn].location[1] - mesh.points[ids].location[1]);
            }
        }

        init.FieldList.insert(std::pair<std::string,Field::field>("gradx"+Scalar,gradx_field));
        init.FieldList.insert(std::pair<std::string,Field::field>("grady"+Scalar,grady_field));


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

//    projection(bool implicit = false){
//        this->implicit = implicit;
//    }

    void predictor(Init &init,Mesh &mesh){
        Field::field Vxstar,Vystar;
        Vxstar.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        Vystar.field.resize(mesh.Y_seed,std::vector<double> (mesh.X_seed));
        this->Conv(init,mesh,"Ux");//modify the location
        this->Conv(init,mesh,"Uy");
        this->Diff(init,mesh,"Ux");
        this->Diff(init,mesh,"Uy");
        //Read file must apply after init action.
        Field::field  &Ux = init.FieldList.find("Ux")->second;
        Field::field  &Uy = init.FieldList.find("Uy")->second;
        Field::field &convUx = init.FieldList.find("convUx")->second;
        Field::field &convUy = init.FieldList.find("convUy")->second;
        Field::field &diffUx = init.FieldList.find("diffUx")->second;
        Field::field &diffUy = init.FieldList.find("diffUy")->second;



        for(int j = 0 ; j < mesh.Y_seed-1;j++){
            for(int i = 0; i < mesh.X_seed-1;i++){
                Vxstar.field[j][i] = Ux.field[j][i] + dt * (-convUx.field[j][i] + mu * diffUx.field[j][i]);
                Vystar.field[j][i] = Uy.field[j][i] + dt * (-convUy.field[j][i] + mu * diffUy.field[j][i]);
            }
        }
        init.FieldList.insert(std::pair<std::string,Field::field>("Vxstar",Vxstar));
        init.FieldList.insert(std::pair<std::string,Field::field>("Vystar",Vystar));

    }

    void Possion(Init &init,Mesh &mesh){

        this->Div(init,mesh,"Vxstar","Vystar");
        std::vector<std::vector<double>> & DivUstar = init.FieldList.find("divU")->second.field;
        std::vector<double> RHS(mesh.Y_seed * mesh.X_seed,0);

        for(int j=0;j<mesh.Y_seed;j++){
            for(int i=0;i<mesh.X_seed;i++)
            {
                RHS[j*mesh.X_seed+i] = DivUstar[j][i] * rho/dt;
            }
        }

        std::vector<std::vector<double>> LaplaceP(mesh.Y_seed*mesh.X_seed,std::vector<double> (mesh.Y_seed*mesh.X_seed,0));
//
        double X_length = std::stof(mesh.FindConfig(mesh.keyValueMap,"X_length"));
        double Y_length = std::stof(mesh.FindConfig(mesh.keyValueMap,"Y_length"));
        double X_seed = std::stof(mesh.FindConfig(mesh.keyValueMap,"X_seed"));
        double Y_seed = std::stof(mesh.FindConfig(mesh.keyValueMap,"Y_seed"));
        double dx = X_length/X_seed;
        double dy = Y_length/Y_seed;

        for(int id = 0; id<mesh.X_seed*mesh.Y_seed;id++){

            std::vector<int> IJK = mesh.IDtoIJK(id);
            int ide = mesh.IJKtoID(IJK[0]+1,IJK[1],0);
            int idw = mesh.IJKtoID(IJK[0]-1,IJK[1],0);
            int idn = mesh.IJKtoID(IJK[0],IJK[1]+1,0);
            int ids = mesh.IJKtoID(IJK[0],IJK[1]-1,0);
//
            if(IJK[0]==0){
                //dpdx=0
                LaplaceP[id][id] += 1/dx;
                LaplaceP[id][ide] -= 1/dx;
            }
            else if(IJK[0] == mesh.X_seed - 1){
                   //pdx =atm
//                LaplaceP[id][id] += 1/dx;
//                LaplaceP[id][idw] -= 1/dx;
                LaplaceP[id][id] = 1;
                double Patm = std::stof(mesh.FindConfig(mesh.keyValueMap,"Pout"));
                RHS[id] = Patm;
            }else if(IJK[1] == 0){
                LaplaceP[id][id] += 1/dy;
                LaplaceP[id][idn] -= 1/dy;
            }else if(IJK[1] == mesh.Y_seed - 1){
                LaplaceP[id][id] += 1/dy;
                LaplaceP[id][ids] -= 1/dy;
            }else{
                /*
                 *                    Pn
                 *         Pw    |    P    |   Pe
                 *                    Ps
                 */
                LaplaceP[id][id] += -2/pow(dx,2)  -2/pow(dy,2);
                LaplaceP[id][idw] += 1/pow(dx,2);
                LaplaceP[id][ide] += 1/pow(dy,2);
                LaplaceP[id][idn] += 1/pow(dy,2);
                LaplaceP[id][ids] += 1/pow(dy,2);
            }

        }

        arma::mat A(mesh.Y_seed*mesh.X_seed,mesh.X_seed*mesh.Y_seed);
        arma::mat b = arma::conv_to<arma::mat>::from(RHS);

        for(int j=0;j<mesh.Y_seed*mesh.X_seed;j++){
            for(int i=0;i<mesh.Y_seed*mesh.X_seed;i++){
                A(j,i) = LaplaceP[j][i];
            }
        }

        arma::mat x = arma::solve(A,b);

        Field::field &P = init.FieldList.find("P")->second;
        for(int i = 0;i <mesh.X_seed*mesh.Y_seed;i++){
            std::vector<int> ijk = mesh.IDtoIJK(i);
            P.field[ijk[1]][ijk[0]] = x(i);
        }


    }

    void corrector(Init &init,Mesh &mesh)
    {
        this->Grad(init,mesh,"P");
        //Read file must apply after init action.
        Field::field  &dpdx = init.FieldList.find("gradxP")->second;
        Field::field  &dpdy = init.FieldList.find("gradyP")->second;
        Field::field  &Ux   = init.FieldList.find("Ux")->second;
        Field::field  &Uy   = init.FieldList.find("Uy")->second;
        Field::field  &Uxstar   = init.FieldList.find("Vxstar")->second;
        Field::field  &Uystar   = init.FieldList.find("Vystar")->second;
        for(int j=1;j<mesh.Y_seed-1;j++){
            for(int i=1;i<mesh.X_seed-1;i++){
                Ux.field[j][i] = Uxstar.field[j][i] - dt/rho * dpdx.field[j][i];
                Uy.field[j][i] = Uystar.field[j][i] - dt/rho * dpdy.field[j][i];
            }
        }
    }


};

#endif //GSOLVER_SOLVER_H
