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


class Point;

class Mesh {

private:
    std::string file_name;

public:
    std::string key,value;
    std::vector<Point> points;

    std::map<std::string,std::string> ReadConfiguration(std::string const&file_name)
    {
        std::stringstream ss;
        std::ifstream Mesh_Config(file_name);
        std::string line;
        std::map<std::string,std::string> key_value_map;

        while (std::getline(Mesh_Config,line)) {
            //Treat as comment after //
            line = std::regex_replace(line, std::regex("//.*"), " ");
            std::cout<<line<<std::endl;
            ss << line << std::endl;
        }

        while(ss>>key>>value){
            key_value_map.insert(std::pair<std::string,std::string>(key,value));
        };
        return key_value_map;
    }

    void GenerateMesh(std::map<std::string,std::string> Config)
    {
        //Get input
        double X_length = std::stof(this->FindConfig(Config,"X_length"));
        double Y_length = std::stof(this->FindConfig(Config,"Y_length"));
        int X_seed = std::stoi(this->FindConfig("X_seed"));
        int Y_seed = std::stoi(this->FindConfig("Y_seed"));

        //if mesh size equal
        double x_size = X_length/(X_seed - 1);
        double y_size = Y_length/(Y_seed - 1);

        //
        int ID=0;

        for(int i=0;i<Y_seed;i++){
            for(int j=0;j<X_seed;j++){
                Point point;
                point.id = ID;
                point.location.push_back(x_size*j);
                point.location.push_back(y_size*i);
                this->points.push_back(point);
            }
        }


    }

    std::string FindConfig(std::map<std::string,std::string> Config,std::string key){
        if(Config.find(key) != Config.end()){
            return Config.find(key)->second;
        }else{
            std::cout<<"Config "+ key + "Not Found!"<<std::endl;
            return "Nan";
        }
    }

    class Point {

    public:
        int id;
        std::vector<double> location;
        std::vector<int> adj;
        bool isBC=false;
        std::string type;
    };

};




#endif //GSOLVER_MESH_H
