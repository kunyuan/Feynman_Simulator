//
//  parameter_map.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__parser__
#define __Feynman_Simulator__parser__

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "utility/utility.h"
#include "utility/abort.h"

#define SetPara(para, value) (para).Set((#value), (value));
#define GetPara(para, value) (para).Get((#value), (value));

/**
*  Parse configuration file with the format:
*
*  Value1          #key1
*  Value2          #key2
*  ##Some comments
*  Value3          #key3
*  ......          #....
*
*  to use the set member function of the class, you have to overload ToString(T) method for your type T
*/
class SimpleParser {
  public:
    bool ParseFile(const std::string &);
    void SaveToFile(const std::string &, std::string Mode = "a");
    std::string PrettyString();

    void clear();

    template <typename T>
    void Set(std::string key, T value)
    {
        _map[key] = ToString(value);
    }
    template <typename T>
    void Set(std::string key, std::vector<T> value)
    {
        ASSERT_ALLWAYS(value.size() != 0, "vector should has element in it!");
        _map[key] = "[" + ToString(value[0]);
        for (auto iter = ++value.begin(); iter < value.end(); ++iter)
            _map[key] += "," + ToString(*iter);
        _map[key] += "]";
    }

    template <typename T>
    void Get(std::string key, T &value)
    {
        _MakeSureKeyExists(key);
        std::stringstream ss(_map.at(key));
        if (ss.peek() == '[') {
            char c;
            ss >> c;
        }
        ss >> value;
        if (ss.fail())
            ABORT("Fail to read " << key << "!");
    }
    template <typename T>
    void Get(std::string key, std::vector<T> &value, char sep = ',')
    {
        _MakeSureKeyExists(key);
        std::stringstream ss(_map.at(key));
        if (ss.peek() == '[') {
            char c;
            ss >> c;
        }
        while (ss.peek() != EOF && ss.peek() != ']') {
            if (ss.peek() != sep) {
                T elem;
                //                std::string temp;
                //                ss >> temp;
                //                std::cout << temp << std::endl;
                ss >> elem;
                if (ss.fail())
                    ABORT("Fail to read "
                          << "," << key << "!");
                value.push_back(elem);
            }
            else
                ss.get();
        }
    }

    std::pair<std::string, std::string> make_pair(std::string);

  private:
    std::map<std::string, std::string> _map;
    bool _MakeSureKeyExists(std::string name);
    bool _DoesKeyExist(std::string name);
};

#endif /* defined(__Feynman_Simulator__parser__) */
