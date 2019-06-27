//
// Created by konstantin on 16.06.19.
//

#ifndef NEWFEM_MANAGER_H
#define NEWFEM_MANAGER_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <cstring>

#include "./structurs/mesh.h"
#include "./structurs/paths.h"
#include "./structurs/arguments.h"

using json = nlohmann::json;

class Manager {
public:
    static Paths paths;
    static Arguments arguments;
    static Mesh mesh;

    static string name;

    static void create(json inputData, string location);
    static void finish();

    static json inputData;

    static void parseMesh(string pathToFolderOfMesh, float h, string nameMesh);
    static string createFolder(string pathToResult, string nameResult, bool force);
    static void writeLog(string path, string nameResult, json arg);
    static void converteToBinary(string pathToFolder, string pathToPreplot);
    static void archiveFolder(string pathToFolder, bool deleteSource);

    static void recordDataTecplot(vector<double*> column, vector<string> nameColumn, float time);

    static int runCommand(string command) {return system (command.c_str());}
};


#endif //NEWFEM_MANAGER_H
