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
#include "./TwoPhaseFlow.h"

#include "./Manager.h"

using namespace std;
using json = nlohmann::json;

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

    Manager::create(inputData, location);

    //SUPG supg;
    NavierStokes ns;
    //TwoPhaseFlow twf;

    int timeB = (int)(clock() / 1000000.0);

    cout << "TIME ALL = " << timeB - timeA << " sec" << endl;

    return 0;
}