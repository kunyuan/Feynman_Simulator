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
#include "utility/pyglue/type_cast_interface.h"
#include "utility/pyglue/pywrapper.h"

#define SET(para, value) (para).Set((#value), (value));
#define GET(para, value) (para).Get((#value), (value));
class Dictionary : public Python::ITypeCast {
public:
    friend class Dictionary;
    Dictionary();
    Dictionary(const Python::Object& obj)
        : _Dict(obj)
    {
        if (!PyDict_Check(obj.Get()))
            ERRORCODEABORT(ERR_VALUE_INVALID, "Expected PyDict object to create Dictionary!");
    }
    Dictionary(const Dictionary& dict)
        : _Dict(dict.GetObject())
    {
    }
    ~Dictionary()
    {
        _Dict.Destroy();
    }

    virtual Python::Object CastToPy() const;
    virtual bool Convert(Python::Object);
    template <typename T>
    void Set(const std::string& key, const T& value)
    {
        Python::AnyObject object(value);
        PyDict_SetItemString(_Dict.Get(), key.c_str(), object.Get());
    }
    template <typename T>
    bool Get(const std::string& key, T& value)
    {
        Python::AnyObject object = Python::Object(PyDict_GetItemString(_Dict.Get(), key.c_str()), Python::NoRef);
        if (object.Get() == nullptr)
            return false;
        value = object.As<T>();
        return true;
    }
    template <typename T>
    T Get(std::string key)
    {
        T value;
        if (Get(key, value) == false)
            ERRORCODEABORT(ERR_KEY_NOT_FOUND,
                           "Key " << key << " is not found in dictionary!");
        return value;
    }
    void Clear();
    void Print();
    std::string PrettyString() { return _Dict.PrettyString(); }
    void LoadFromString(const std::string&);
    void Load(const std::string& FileName);
    void Save(const std::string& FileName, const std::string& Mode = "a");

    Python::AnyObject GetObject() const { return _Dict; }
    void _PrintDebug() const;

private:
    Python::AnyObject _Dict;
};
int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
