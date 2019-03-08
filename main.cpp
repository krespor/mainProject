#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <cstring>
#include <time.h>

#include "./structurs/mesh.h"
#include "./structurs/paths.h"
#include "./structurs/arguments.h"

#include "./NavierStokes.h"
#include "./SUPG.h"

using namespace std;
using json = nlohmann::json;

Mesh parseMesh(string pathToFolderOfMesh, float h, string nameMesh);
string  createFolder(string pathToResult, string nameResult);
void converteToBinary(string pathToFolder, string pathToPreplot);
void archiveFolder(string pathToFolder, bool deleteSource);
void createLogFolder(string pathToLog, string nameFolder);
void writeLog(string path, string nameMesh);

int runCommand(string command) {return system (command.c_str());}

int main(int argc, char *argv[])
{
    int timeA = (int)(clock() / 1000000.0);

    string pathToInput = "/home/konstantin/Solver/inputData.json";
    string location;

    if (argc > 1)
    {
        pathToInput = argv[1];
        location = "claster";
    } else
    {
        location = "local";
    }

    ifstream ifs(pathToInput);
    json inputData = json::parse(ifs);
    ifs.close();

    Paths paths;
    Arguments arguments;

    paths.mesh = inputData["paths"][location]["mesh"];
    paths.result = inputData["paths"][location]["result"];
    paths.log = inputData["paths"][location]["log"];
    paths.preplot = inputData["paths"][location]["preplot"];

    arguments.write.interval = inputData["arguments"]["write"]["interval"];
    arguments.write.onlyLastFrame = inputData["arguments"]["write"]["saveOnlyLastFrame"];
    arguments.write.record = inputData["arguments"]["write"]["record"];

    arguments.time.del_t = inputData["arguments"]["time"]["del_t"];
    arguments.time.start = inputData["arguments"]["time"]["start"];
    arguments.time.end = inputData["arguments"]["time"]["end"];

    arguments.individualArguments = inputData["arguments"]["individual"];

    Mesh mesh = parseMesh(paths.mesh, inputData["h"], inputData["mesh"]);

    if (arguments.write.record)
        mesh.name = createFolder(paths.result, inputData["name"]);
    paths.result += mesh.name + "/";

    createLogFolder(paths.log, mesh.name);

    writeLog(paths.log, mesh.name);

    //SUPG *supg = new SUPG(mesh, arguments, paths.result);
    //delete supg;

    NavierStokes *ns = new NavierStokes(mesh, arguments, paths.result);
    delete ns;

    if (inputData["replaceToBinary"])
        if(arguments.write.record)
            converteToBinary(paths.result, paths.preplot);

    if (inputData["archive"])
        if(arguments.write.record)
            archiveFolder(paths.result, inputData["deleteFolderAfterArchive"]);

    int timeB = (int)(clock() / 1000000.0);

    cout << "ttt = " << timeB - timeA << endl;

    return 0;
}

void archiveFolder(string pathToFolder, bool deleteSource)
{
    pathToFolder.erase(strlen(pathToFolder.c_str()) - 1);

    if (runCommand("tar -cf " + pathToFolder + ".tar.gz " + pathToFolder) == 0 && deleteSource)
    {
        runCommand("rm -rf " + pathToFolder);
    }

}

void writeLog(string path, string nameMesh)
{
    json log =
            {
                    {"nameMesh", nameMesh}
            };
    ofstream file(path + nameMesh + "/" + "log.json");
    file << log << endl;
    file.close();
}

void converteToBinary(string pathToFolder, string pathToPreplot)
{
    string currentDatFile, prevDatFile;
    string nameListFiles = "listFile.txt";

    ifstream file;

    runCommand("ls " + pathToFolder + " > " + nameListFiles);

    file.open(nameListFiles);

    prevDatFile = "";
    file >> currentDatFile;

    while (currentDatFile != prevDatFile)
    {
        prevDatFile = currentDatFile;
        if (runCommand(pathToPreplot + " " + pathToFolder + currentDatFile) == 0)
        {
            runCommand("rm -f " + pathToFolder + currentDatFile);
        }

        file >> currentDatFile;
    }
    runCommand("rm -rf " + nameListFiles);
    file.close();
}

string createFolder(string pathToResult, string nameResult)
{
    string tmp = pathToResult + nameResult;
    string tmpUnix = "mkdir " + tmp;
    int i;
    for (i = 0; true; i++)
    {
        if (system(tmpUnix.c_str()) != 0)
        {
            tmpUnix = "mkdir " + tmp;
            tmpUnix += "[" + to_string(i) + "]";
        } else
        {
            if (i == 0)
                pathToResult += nameResult + "/";
            else
               pathToResult += nameResult + "[" + to_string(i - 1) + "]/";
            break;
        }
    }
    if (i != 0)
    {
        return nameResult + "[" + to_string(i - 1) + "]";
    }
    return nameResult;
}

void createLogFolder(string pathToLog, string nameFolder)
{
    pathToLog += nameFolder;
    string tmpUnix = "mkdir " + pathToLog;

    if (system(tmpUnix.c_str()) != 0)
    {
        tmpUnix = "rm -rf " + pathToLog;
        system(tmpUnix.c_str());
        tmpUnix = "mkdir " + pathToLog;
        system(tmpUnix.c_str());
    }
}

Mesh parseMesh(string pathToFolderOfMesh, float h, string nameMesh)
{
    Mesh mesh;

    mesh.h = h;

    ostringstream strStream;
    strStream << h;
    string hStr(strStream.str());
    hStr[1] = '_';

    string pathToMesh = pathToFolderOfMesh + nameMesh + "/" + hStr + "/" + nameMesh;

    ifstream file;


    file.open(pathToMesh + ".dat");
    if (!file)
    {
        throw invalid_argument( "ERROR: mesh '" + pathToMesh + ".dat' dont exist!");
    }

    double musor;
    unsigned int *arrayMusor = new unsigned int[4];
    int tmpM = 0;

    file >> mesh.n >> tmpM;
    mesh.nodes = new double*[mesh.n];
    for (unsigned int i = 0; i < mesh.n; i++)
    {
        mesh.nodes[i] = new double[2];
    }

    for (unsigned int i = 0; i < mesh.n; i++)
    {
        file >> musor >> mesh.nodes[i][0] >> mesh.nodes[i][1] >> musor;
    }
    for (unsigned int i = 0; true; i++)
    {

        file >> arrayMusor[0] >> arrayMusor[1] >> arrayMusor[2] >> arrayMusor[3];
        if (arrayMusor[1] != 102)
        {
            mesh.m = tmpM - i;
            break;
        }
    }
    mesh.elements = new unsigned int*[mesh.m];
    for (unsigned int i = 0; i < mesh.m; i++)
    {
        mesh.elements[i] = new unsigned int[3];
    }

    mesh.elements[0][0] = arrayMusor[2];
    mesh.elements[0][1] = arrayMusor[3];
    file >> mesh.elements[0][2];
    for (unsigned int i = 1; i < mesh.m; i++)
    {
        file >> musor >> musor >> mesh.elements[i][0] >> mesh.elements[i][1] >> mesh.elements[i][2];
    }
    for (unsigned int i = 0; i < mesh.m; i++)
    {
        mesh.elements[i][0]--;
        mesh.elements[i][1]--;
        mesh.elements[i][2]--;
    }
    delete[]arrayMusor;
    //////////////////////////////////////////////////////////////////////////
    file.close();

    mesh.countBorder = 0;

    for (char i = '0'; true; i++)
    {
        file.open(pathToMesh + "G" + i + ".dat");
        if (!file)
        {
            if (i == '0')
            {
                throw invalid_argument( "ERROR: border '" + pathToMesh + "G" + i + ".dat" + ".dat' dont exist!");
            }
            break;
        }
        mesh.countBorder++;
        file.close();
    }

    mesh.border = new Border[mesh.countBorder];

    for (int i = 0; i < mesh.countBorder; i++)
    {
        file. open(pathToMesh + "G" + to_string(i) + ".dat");

        file >> mesh.border[i].n >> mesh.border[i].m;

        mesh.border[i].numberNodes = new unsigned int[mesh.border[i].n];
        mesh.border[i].line = new unsigned int*[mesh.border[i].m];
        for (unsigned int j = 0; j < mesh.border[i].m; j++)
        {
            mesh.border[i].line[j] = new unsigned int[2];
        }

        double musor;
        for (unsigned int j = 0; j < mesh.border[i].n; j++)
        {
            file >> mesh.border[i].numberNodes[j] >> musor >> musor >> musor;
            mesh.border[i].numberNodes[j]--;
        }
        for (unsigned int j = 0; j < mesh.border[i].m; j++)
        {
            file >> musor >> musor >> mesh.border[i].line[j][0] >> mesh.border[i].line[j][1];
            mesh.border[i].line[j][0]--;
            mesh.border[i].line[j][1]--;
        }

        file.close();
    }

    return mesh;
}