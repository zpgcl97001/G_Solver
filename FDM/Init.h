//
// Created by Guochenlin on 12/24/2021.
//

#ifndef GSOLVER_INIT_H
#define GSOLVER_INIT_H

#include "Mesh.h"
#include "Field.h"
#include <math.h>

class Init {

private:

public:

    std::map<std::string,Field::field> FieldList;

    void MammalInit(Mesh &mesh) {
        //Init Boundary
        //X = 0  Possunique inflow
        std::vector<Mesh::Point> inlet, outlet, wall, Bluntbody;

        for (Mesh::Point point : mesh.points) {
            if (point.I == 0) {
                inlet.push_back(point);
                point.isBC = true;
                point.type = "inlet";
            } else if (point.I == mesh.X_seed - 1) {
                outlet.push_back(point);
                point.isBC = true;
                point.type = "outlet";
            } else if (point.J == 0 || point.J == mesh.Y_seed - 1) {
                wall.push_back(point);
                point.isBC = true;
                point.type = "wall";
            } else if (mesh.isConfig(mesh.keyValueMap, "Bluntbody")) {
                double bl_x_lo = std::stof(mesh.FindConfig(mesh.keyValueMap, "Blx-"));
                double bl_y_lo = std::stof(mesh.FindConfig(mesh.keyValueMap, "Bly-"));
                double bl_x_up = std::stof(mesh.FindConfig(mesh.keyValueMap, "Blx+"));
                double bl_y_up = std::stof(mesh.FindConfig(mesh.keyValueMap, "Bly+"));
                if (point.location[0] >= bl_x_lo && point.location[0] <= bl_x_up &&
                    point.location[1] >= bl_y_lo && point.location[1] <= bl_y_up) {
                    point.isBC = true;
                    point.type = "wall";
                    Bluntbody.push_back(point);
                }
            }

        }

        mesh.BoundaryList.insert(std::pair<std::string, std::vector<Mesh::Point>>("inlet", inlet));
        mesh.BoundaryList.insert(std::pair<std::string, std::vector<Mesh::Point>>("outlet", outlet));
        mesh.BoundaryList.insert(std::pair<std::string, std::vector<Mesh::Point>>("wall", wall));
        mesh.BoundaryList.insert(std::pair<std::string, std::vector<Mesh::Point>>("Bluntbody", Bluntbody));
        std::cout << "Mesh type init fin." << std::endl;

        //init field
        //Ux
        this->FieldList.insert(std::pair<std::string,Field::field>(this->initField("Ux",mesh)));
        this->FieldList.insert(std::pair<std::string,Field::field>(this->initField("Uy",mesh)));
        this->FieldList.insert(std::pair<std::string,Field::field>(this->initField("Uz",mesh)));
        this->FieldList.insert(std::pair<std::string,Field::field>(this->initField("P",mesh)));

        //Give Boundary condition
        //inlet
        Field::field & Ux = FieldList.find("Ux")->second;
        Field::field & Uy = FieldList.find("Uy")->second;
        Field::field & Uz = FieldList.find("Uz")->second;
        Field::field & P  = FieldList.find("P")->second;
        for(auto point : mesh.BoundaryList.find("inlet")->second){
            //U
            double Um = std::stof(mesh.FindConfig(mesh.keyValueMap,"Ux"));
            double Y_length = std::stof(mesh.FindConfig(mesh.keyValueMap,"Y_length"));
            Ux.field[point.J][point.I] =  Um *(1 - std::pow((double(2*point.J-mesh.Y_seed))/double(mesh.Y_seed),2));
            //P simplify zero Gradient
                //pressure propogate from outlet
            P.field[point.J][point.I] = P.field[point.J][point.I+1];
        }
        for(auto point : mesh.BoundaryList.find("outlet")->second){
            //U extrapolation
            Ux.field[point.J][point.I] = Ux.field[point.J][point.I-1];
            //P = Pout
            double Pout = std::stof(mesh.FindConfig(mesh.keyValueMap,"Pout"));
            P.field[point.J][point.I] = Pout;
        }
        for(auto point : mesh.BoundaryList.find("wall")->second){
            //U no-slip
            Ux.field[point.J][point.I] = 0;
            Uy.field[point.J][point.I] = 0;
            Uz.field[point.J][point.I] = 0;
            P.field[point.J][point.I] = point.J==0? P.field[point.J+1][point.I]:P.field[point.J-1][point.I];
        }

    }



    std::pair<std::string,Field::field> initField(std::string field_name,Mesh &mesh)
    {
        double init_value = std::stof(mesh.FindConfig(mesh.keyValueMap,field_name));
        Field::field field;
        field.name = field_name;
        for(int i=0;i<mesh.Y_seed;i++){
            field.field.emplace_back(mesh.X_seed,0);
            for(int j=0;j<mesh.X_seed;j++){
                field.field[i][j] = init_value;
            }
        }
        return std::pair<std::string,Field::field>(field_name,field);
    }


    };



#endif //GSOLVER_INIT_H
