//
//  parameter_map.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__parser__
#define __Feynman_Simulator__parser__

#include <map>
#include <string>
#include "../utility/utility.h"
#include "../utility/abort.h"

#define SetPara(para, value) (para).set((#value), (value));
#define SetParaArray(para, value, num) (para).set((#value), (value), (num));
#define GetPara(para, value) (para).get((#value), (value));
#define GetParaArray(para, value, num) (para).get((#value), (value), (num));

/**
*  to use the set member function of the class, you have to overload ToString(T) method for your type T
*/
class SimpleParser {
  public:
    bool ParseFile(const std::string &, bool AbortIfFail = true);
    void SaveToFile(const std::string &, std::string Mode = "a");
    std::string PrettyString();

    void clear();

    void addKey(std::string key);
    void eraseKey(std::string key);

    template <typename T>
    void set(std::string key, T value)
    {
        key = _ToUpper(key);
        _map[key] = ToString(value);
    }

    template <typename T>
    void set(std::string key, T *value, int num)
    {
        key = _ToUpper(key);
        _map[key] = ToString(value[0]);
        for (int i = 1; i < num; i++)
            _map[key] += ("," + ToString(value[i]));
    }

    std::string get(std::string key);
    template <typename T>
    void get(std::string key, T &value)
    {
        key = _ToUpper(key);
        _MakeSureKeyExists(key);
        std::stringstream ss(_map.at(key));
        ss >> value;
        if (ss.fail())
            ABORT("Fail to read " << key << "!");
    }

    template <typename T>
    void get(std::string key, T *value, int num, char sep = ',')
    {
        key = _ToUpper(key);
        _MakeSureKeyExists(key);
        std::stringstream ss(_map.at(key));

        char sepchar;
        ss >> value[0];
        for (int i = 1; i < num; i++) {
            ss >> sepchar;
            if (sepchar != sep)
                ABORT("Sep char " << sepchar << " is not expected as the separator. I will expect " << sep);
            ss >> value[i];
        }
        if (ss.fail())
            ABORT("Fail to read " << key << "!");
    }

    std::pair<std::string, std::string> make_pair(std::string);

  private:
    std::map<std::string, std::string> _map;
    std::string _ToUpper(std::string source);
    bool _MakeSureKeyExists(std::string name);
    bool _DoesKeyExist(std::string name);
};

#endif /* defined(__Feynman_Simulator__parser__) */
