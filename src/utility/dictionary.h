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
#include "utility/abort.h"
#include "utility/pyglue/pywrapper.h"

#define SET(para, value) (para[#value] = (value));
#define GET(para, value) (para).Get((#value), (value));
typedef std::map<std::string, Python::AnyObject> PythonMap;
class Dictionary : public Python::ITypeCast {
private:
    PythonMap _Map;

public:
    friend class Dictionary;
    Dictionary()
    {
    }
    Dictionary(const std::string& script)
    {
        LoadFromString(script);
    }
    Dictionary& operator=(const Dictionary& dict)
    {
        _Map = dict._Map;
        return *this;
    }

    //ITypeCast interface
    virtual Python::Object ToPy() const;
    virtual bool FromPy(const Python::Object&);
    Python::AnyObject& operator[](const std::string& key);

    template <typename T>
    Dictionary& Set(const std::string& key, const T& value)
    {
        _Map[key] = Python::AnyObject(value);
        return *this;
    }
    template <typename T>
    bool Get(const std::string& key, T& value) const
    {
        return Python::Convert(_Map.at(key), value);
    }
    template <typename T>
    T Get(const std::string& key) const
    {
        T value;
        if (!Python::Convert(_Map.at(key), value))
            ERRORCODEABORT(ERR_VALUE_INVALID, "Fail to convert " << key);
        return value;
    }
    bool HasKey(const std::string& key) const;
    void Clear();
    bool IsEmpty() const;
    void LoadFromString(const std::string&);
    void Load(const std::string& FileName, const std::string& key = "default");
    void Save(const std::string& FileName, const std::string& Mode = "a",
              const std::string& key = "default");
    void BigLoad(const std::string& FileName);
    void BigSave(const std::string& FileName);
    void Print() const;
    std::string PrettyString() const;

    const PythonMap::iterator begin()
    {
        return _Map.begin();
    }
    const PythonMap::iterator end()
    {
        return _Map.end();
    }
};
int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
