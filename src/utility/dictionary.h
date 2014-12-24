//
//  serialization.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__serialization__
#define __Feynman_Simulator__serialization__

//#include "numpy/arrayobject.h"

#include "utility/convention.h"
#include "utility/utility.h"
#include "utility/vector.h"
#include "utility/abort.h"
#include <map>
#include "utility/pyglue/pywrapper.h"

#define SET(para, value) (para).Set((#value), (value));
#define GET(para, value) (para).Get((#value), (value));
class Dictionary : public Python::ITypeCast {
private:
    std::map<std::string, Python::AnyObject> _Map;

public:
    friend class Dictionary;
    Dictionary()
    {
    }
    Dictionary& operator=(const Dictionary& dict)
    {
        _Map = dict._Map;
        return *this;
    }

    //ITypeCast interface
    virtual Python::Object ToPy() const;
    virtual bool FromPy(Python::Object);

    template <typename T>
    void Set(const std::string& key, const T& value)
    {
        _Map[key] = value;
    }
    template <typename T>
    bool Get(std::string key, T& value)
    {
        return Python::Convert(_Map[key], value);
    }
    template <typename T>
    T Get(std::string key)
    {
        return _Map[key].As<T>();
    }
    void Clear();
    void LoadFromString(const std::string&);
    void Load(const std::string& FileName);
    //the key will be used as the name of Dictionary when Mode="a"
    void Save(const std::string& FileName, const std::string& Mode = "a");
    void Print();
    std::string PrettyString();
};
int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
