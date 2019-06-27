//
// Created by konstantin on 16.06.19.
//

#include "Manager.h"

Paths Manager::paths;
Arguments Manager::arguments;
Mesh Manager::mesh;

json Manager::inputData;

string Manager::name;

void Manager::create(json inputData, string location)
{
    Manager::inputData = inputData;

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

    parseMesh(paths.mesh, inputData["h"], inputData["mesh"]);

    if (arguments.write.record)
    {
        name = createFolder(paths.result, inputData["name"], inputData["forceMode"]);
        writeLog(paths.log, name, inputData);
    }

    paths.result += name + "/";
}

void Manager::parseMesh(string pathToFolderOfMesh, float h, string nameMesh)
{
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
}

string Manager::createFolder(string pathToResult, string nameResult, bool force) {
    string tmp = pathToResult + nameResult;
    string tmpUnix = "mkdir " + tmp;
    int i;
    for (i = 0; true; i++)
    {
        if (system(tmpUnix.c_str()) != 0)
        {
            if (force)
            {
                tmpUnix = "rm -rf " + tmp;
                system(tmpUnix.c_str());
                tmpUnix = "mkdir " + tmp;
            } else
            {
                tmpUnix = "mkdir " + tmp;
                tmpUnix += "[" + to_string(i) + "]";
            }
        } else
        {
            if (i == 0 || force)
                pathToResult += nameResult + "/";
            else
                pathToResult += nameResult + "[" + to_string(i - 1) + "]/";
            break;
        }
    }
    if (i == 0 || force)
    {
        return nameResult;
    }
    return nameResult + "[" + to_string(i - 1) + "]";
}

void Manager::writeLog(string path, string nameResult, json arg)
{
    json log =
            {
                    {"h", arg["h"]},
                    {"name", nameResult},
                    {"mesh", arg["mesh"]},
                    {"interval", arg["arguments"]["write"]["interval"]}
            };
    ofstream file(path + nameResult + ".json");
    file << log << endl;
    file.close();
}

void Manager::finish()
{
    if (inputData["replaceToBinary"])
        if(arguments.write.record)
            converteToBinary(paths.result, paths.preplot);

    if (inputData["archive"])
        if(arguments.write.record)
            archiveFolder(paths.result, inputData["deleteFolderAfterArchive"]);
}

void Manager::converteToBinary(string pathToFolder, string pathToPreplot)
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
        if (runCommand(pathToPreplot + " " + pathToFolder + currentDatFile + ">> /dev/null") == 0)
        {
            runCommand("rm -f " + pathToFolder + currentDatFile);
        }

        file >> currentDatFile;
    }
    runCommand("rm -rf " + nameListFiles);
    file.close();
}

void Manager::archiveFolder(string pathToFolder, bool deleteSource)
{
    pathToFolder.erase(strlen(pathToFolder.c_str()) - 1);

    if (runCommand("tar -cf " + pathToFolder + ".tar " + pathToFolder) == 0 && runCommand("gzip  " + pathToFolder + ".tar ") == 0 && deleteSource)
    {
        runCommand("rm -rf " + pathToFolder);
    }
}

void Manager::recordDataTecplot(vector<double *> column, vector<string> nameColumn, float time)
{
    ofstream file;

    file.open(paths.result + to_string(time) + ".dat");

    file << "VARIABLES = X, Y,";
    file.precision(12);

    for (string nc : nameColumn)
    {
        file << " ," << nc;
    }

    file << endl;
    file << "ZONE T=\"" + to_string(time)  + "\", N=" << mesh.n << ", E=" << mesh.m << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE" << endl;

    for (unsigned int i = 0; i < mesh.n; i++)
    {
        file << mesh.nodes[i][0] << " " << mesh.nodes[i][1] << " ";
        for (unsigned int j = 0; j < column.size(); j++)
        {
            file << column[j][i] << " ";
        }
        file << endl;
    }

    for (unsigned int i = 0; i < mesh.m; i++)
    {
        file << mesh.elements[i][0] + 1 << " " << mesh.elements[i][1] + 1 << " " << mesh.elements[i][2] + 1 << endl;
    }
    file.close();
}