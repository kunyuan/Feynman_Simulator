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
#include <iostream>

using std::string;

namespace Python {

bool Convert(Object obj, AnyObject& any)
{
    any = obj;
    return true;
}
void AnyObject::EvalScript(const std::string& script)
{
    PyGILState_STATE state = PyGILState_Ensure();
    ON_SCOPE_EXIT([&] {PyGILState_Release(state); });
    Object main = Object(PyImport_AddModule("__main__"), NoRef);
    Object global = Object(PyModule_GetDict(main.Get()), NoRef); //borrowed reference
    Object local = PyDict_New();
    *this = PyRun_String(script.c_str(), Py_eval_input,
                         global.Get(), local.Get());
    PropagatePyError();
}

ModuleObject::ModuleObject(const Object& obj)
    : Object(obj)
{
    if (!PyModule_Check(obj.Get()))
        ABORT("PyModule object is expected!");
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
    PropagatePyError();
    ASSERT_ALLWAYS(module.Get() != nullptr, "Failed to import module " << script_path);
    *this = module;
}

Object ModuleObject::CallFunction(const std::string& name)
{
    PyGILState_STATE state = PyGILState_Ensure();
    ON_SCOPE_EXIT([&] {PyGILState_Release(state); });
    Object func = GetAttr(name);
    ASSERT_ALLWAYS(PyCallable_Check(func.Get()), name << " is not a callable function!");
    Object ret = PyObject_CallObject(func.Get(), nullptr);
    PropagatePyError();
    return { ret };
}

Object ModuleObject::GetAttr(const std::string& name)
{
    Object obj = PyObject_GetAttrString(_PyPtr, name.c_str());
    PropagatePyError();
    ASSERT_ALLWAYS(obj.Get() != nullptr, "Failed to get attribute: " << name);
    return { obj };
}

bool ModuleObject::HasAttr(const std::string& name)
{
    try {
        GetAttr(name);
        return true;
    }
    catch (KeyInvalid) {
        return false;
    }
}
}