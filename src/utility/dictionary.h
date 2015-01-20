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

#define SET(para, value) (para[#value] = (value))
#define GET(para, value) para.Get((#value), (value))
#define GET_WITH_DEFAULT(para, value, default) \
    do {                                       \
        if (para.HasKey(#value))               \
            para.Get((#value), (value));       \
        else                                   \
            value = default;                   \
    } while (0)
typedef std::map<std::string, Python::AnyObject> PythonMap;
class Dictionary : public Python::ITypeCast {
private:
    PythonMap _Map;

public:
    Dictionary()
    {
    }
    Dictionary(const std::string& script)
    {
        LoadFromString(script);
    }
    template <typename T>
    Dictionary(const std::string& key, const T& value)
    {
        _Map[key] = value;
    }
    Dictionary& operator=(const Dictionary& dict)
    {
        if (this == &dict)
            return *this;
        _Map = dict._Map;
        return *this;
    }

    Python::AnyObject& operator[](const std::string& key);
    template <typename T>
    T Get(const std::string& key) const
    {
        T value;
        if (!HasKey(key))
            THROW_ERROR(KeyInvalid, key + " key does not exist!");
        if (!Python::Convert(_Map.at(key), value))
            ABORT("Fail to convert " << key);
        return value;
    }
    template <typename T>
    void Get(const std::string& key, T& value) const
    {
        value = Get<T>(key);
    }

    void Update(const Dictionary&);
    bool HasKey(const std::string& key) const;
    void Clear();
    bool IsEmpty() const;
    void LoadFromString(const std::string&);
    void Load(const std::string& FileName);
    void Save(const std::string& FileName, const std::string& Mode = "a");
    void BigLoad(const std::string& FileName);
    void BigSave(const std::string& FileName);
    void Print() const;
    std::string PrettyString() const;

    //range based iteration
    PythonMap::iterator begin()
    {
        return _Map.begin();
    }
    PythonMap::const_iterator begin() const
    {
        return _Map.begin();
    }
    PythonMap::iterator end()
    {
        return _Map.end();
    }
    PythonMap::const_iterator end() const
    {
        return _Map.end();
    }
    //ITypeCast interface
    virtual Python::Object ToPy() const;
    virtual bool FromPy(const Python::Object&);
};
int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
