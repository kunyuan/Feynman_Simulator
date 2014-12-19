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
#include "utility/pyglue/pywrapper.h"

class Dictionary {
public:
    friend class Dictionary;
    Dictionary();
    //borrow reference, no change on reference of PyObject
    Dictionary(const Python::Object& obj)
        : _Dict(obj)
    {
    }
    //borrow reference, no change on reference of PyObject
    Dictionary(const Dictionary& dict)
    {
        //        dict._PrintDebug();
        *this = dict;
        //        dict._PrintDebug();
    }
    template <typename T>
    void Set(const std::string& key, const T& value)
    {
        Python::Object object(value);
        PyDict_SetItemString(_Dict.get(), key.c_str(), object.get());
    }
    void Set(const std::string& key, const Dictionary& value)
    {
        PyDict_SetItemString(_Dict.get(), key.c_str(), value._Dict.get());
    }
    template <typename T>
    bool Get(const std::string& key, T& value)
    {
        Python::Object object = PyDict_GetItemString(_Dict.get(), key.c_str());
        if (object.get() == nullptr)
            return false;
        object.Convert(value);
        return true;
    }
    bool Get(const std::string& key, Dictionary& value)
    {
        PyObject* object = PyDict_GetItemString(_Dict.get(), key.c_str());
        if (object == nullptr)
            return false;
        value = Dictionary(object);
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
    void Print();
    void LoadByEval(const std::string&);
    void Load(const std::string& FileName);
    void Save(const std::string& FileName, const std::string& Mode = "a");

private:
    Python::Object _Dict;
    void _PrintDebug() const;
};

int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
