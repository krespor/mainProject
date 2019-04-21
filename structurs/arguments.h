//
// Created by konstantin on 26.02.19.
//

#ifndef NEWFEM_PARAMETRS_H
#define NEWFEM_PARAMETRS_H

#include <nlohmann/json.hpp>

struct Time
{
    float start;
    float end;
    float del_t;
};

struct Write
{
    bool record;
    float interval;
    bool onlyLastFrame;
};

struct Arguments
{
    Time time;
    Write write;
    nlohmann::json individualArguments;
};



#endif //NEWFEM_PARAMETRS_H
