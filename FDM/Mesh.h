//
// Created by Guochenlin on 12/22/2021.
//

#ifndef GSOLVER_MESH_H
#define GSOLVER_MESH_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <vector>
#include <map>
//#include "Mesh.h"


class Mesh {

private:
    std::string file_name;

public:

    class Point {

    public:
        int id;
        std::vector<double> location;
        std::vector<int> adj;
        bool isBC=false;
        std::string type;
        int I,J,K;
    };


    std::string key,value;
    std::vector<Point> points;
    std::map<std::string,std::string> keyValueMap;
    int X_seed,Y_seed;
    std::map<std::string,std::vector<Point>> BoundaryList;



    std::map<std::string,std::string> ReadConfiguration(std::string const&file_name)
    {
        std::stringstream ss;
        std::ifstream Mesh_Config(file_name);
        std::string line;
        std::map<std::string,std::string> key_value_map;

        while (std::getline(Mesh_Config,line)) {
            //Treat as comment after //
            line = std::regex_replace(line, std::regex("//.*"), " ");
            ss << line << std::endl;
        }

        while(ss>>key>>value){
            key_value_map.insert(std::pair<std::string,std::string>(key,value));
        };
        this->keyValueMap = key_value_map;
        return key_value_map;
    }

    void GenerateMesh(std::map<std::string,std::string> Config) {
        //Get input
        double X_length = std::stof(this->FindConfig(Config, "X_length"));
        double Y_length = std::stof(this->FindConfig(Config, "Y_length"));

        this->X_seed = std::stoi(this->FindConfig(Config,"X_seed"));
        this->Y_seed = std::stoi(this->FindConfig(Config,"Y_seed"));



        //if mesh size equal
        double x_size = X_length/(X_seed - 1);
        double y_size = Y_length/(Y_seed - 1);

        //
        int ID=0;

        for(int i=0;i<Y_seed;i++){
            for(int j=0;j<X_seed;j++){
                Point point;
                point.id = ID;
                point.I = j;
                point.J = i;
                point.K = 0;
                point.location.push_back(x_size*j);
                point.location.push_back(y_size*i);
                this->points.push_back(point);
                ID++;
            }
        }


        //ReadBoundary:wall,inlet,outlet,only

    }

    std::string FindConfig(std::map<std::string,std::string> Config,std::string key){
        if(Config.find(key) != Config.end()){
            return Config.find(key)->second;
        }else{
            std::cout<<"Config "+ key + " Not Found!"<<std::endl;
            return "Nan";
        }
    }

    bool isConfig(std::map<std::string,std::string> Config,std::string key){
        if(Config.find(key) != Config.end()){
            return true;
        }else{
            return false;
        }
    }

    std::vector<int> IDtoIJK(int ID){
        int I = ID%this->X_seed;
        int J = ID/this->Y_seed;
        int K = 0;
        std::vector<int> IJK = {I,J,K};
        return IJK;
    }
    int IJKtoID(int I,int J,int K=0){
        return this->Y_seed*this->X_seed*K + this->X_seed*J + I;
    }


};




#endif //GSOLVER_MESH_H
