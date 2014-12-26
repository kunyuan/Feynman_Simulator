/*      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 3 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 *      
 *      Author: 
 *      Matias Fontanini
 * 
 */

#include "utility/utility.h"
#include "utility/complex.h"
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include "pywrapper.h"

using std::string;

namespace Python {

bool Convert(Object obj, AnyObject& any)
{
    any = obj;
    return true;
}
void AnyObject::EvalScript(const std::string& script)
{
    Object main = Object(PyImport_AddModule("__main__"), NoRef);
    Object global = Object(PyModule_GetDict(main.Get()), NoRef); //borrowed reference
    Object local = PyDict_New();
    *this = PyRun_String(script.c_str(), Py_eval_input,
                         global.Get(), local.Get());
    MakeSureNoPyError(ERR_GENERAL);
}

ModuleObject::ModuleObject(const Object& obj)
    : Object(obj)
{
    if (!PyModule_Check(obj.Get()))
        ERRORCODEABORT(ERR_VALUE_INVALID, "PyModule object is expected!");
}
ModuleObject& ModuleObject::operator=(const ModuleObject& obj)
{
    Object::operator=(obj);
    return *this;
}
void ModuleObject::LoadModule(const string& script_path)
{
    char arr[] = "path";
    Object path = Object(PySys_GetObject(arr), NoRef);
    string base_path("."), file_path;
    size_t last_slash(script_path.rfind("/"));
    if (last_slash != string::npos) {
        if (last_slash >= script_path.size() - 2)
            ABORT("Invalid script path");
        base_path = script_path.substr(0, last_slash);
        file_path = script_path.substr(last_slash + 1);
    }
    else
        file_path = script_path;
    if (file_path.rfind(".py") == file_path.size() - 3)
        file_path = file_path.substr(0, file_path.size() - 3);
    Object pwd = PyString_FromString(base_path.c_str());

    PyList_Append(path.Get(), pwd.Get());
    /* We don't need that string value anymore, so deref it */
    Object module = PyImport_ImportModule(file_path.c_str());
    MakeSureNoPyError(ERR_GENERAL);
    *this = module;
}

Object ModuleObject::load_function(const std::string& name)
{
    Object obj = PyObject_GetAttrString(_PyPtr, name.c_str());
    obj.MakeSureNotNull();
    MakeSureNoPyError(ERR_GENERAL);
    return obj;
}

Object ModuleObject::CallFunction(const std::string& name)
{
    Object func = load_function(name);
    //    Object func = load_function(name);
    Object ret = PyObject_CallObject(func.Get(), 0);
    MakeSureNoPyError(ERR_GENERAL);
    return { ret };
}

Object ModuleObject::GetAttr(const std::string& name)
{
    Object obj = PyObject_GetAttrString(_PyPtr, name.c_str());
    MakeSureNoPyError(ERR_GENERAL);
    return { obj };
}

bool ModuleObject::HasAttr(const std::string& name)
{
    try {
        GetAttr(name);
        return true;
    }
    catch (ERRORCODE e) {
        return false;
    }
}
}