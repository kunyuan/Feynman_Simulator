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

#ifndef PYWRAPPER_H
#define PYWRAPPER_H

#include "utility/abort.h"
#include "type_cast.h"
#include "pyarraywrapper.h"
#include "object.h"
#include "utility/scopeguard.h"

namespace Python {
class AnyObject;
bool Convert(Object, AnyObject&);
class AnyObject : public Object {
public:
    AnyObject()
        : Object()
    {
    }
    AnyObject(const Object& obj)
        : Object(obj)
    {
    }
    AnyObject(PyObject* obj, OwnerShip ownership = NewRef)
        : Object(obj, ownership)
    {
    }
    AnyObject& operator=(const AnyObject& obj)
    {
        Object::operator=(obj);
        return *this;
    }

    template <typename T>
    AnyObject(T value)
        : Object(CastToPy(value))
    {
    }
    template <typename T>
    T As()
    {
        T value;
        if (!Python::Convert(*this, value))
            ABORT("Fail to convert PyObject!");
        return value;
    }
    /**
         * \brief Constructs a Python::Object from a script string.
         * 
         * *this be the evaluation of the
         * script. If any errors are encountered while loading this 
         * script, an exception is thrown.
         * 
         * \param script The string of the script to be evaluated.
         */
    void EvalScript(const std::string& script);
};

/**
    *  After ModuleObject construction, you need to call LoadModule to fill it
    */

class ModuleObject : Object {
public:
    /**
         * \brief Constructs a default python object
         */
    ModuleObject()
        : Object()
    {
    }
    ModuleObject(const Object& obj);
    ModuleObject& operator=(const ModuleObject& obj);

    /**
         * \brief Constructs a Python::Object from a script file.
         * 
         * *this will be the representation of the loaded
         * script. If any errors are encountered while loading this 
         * script, an ERRORCODE is thrown.
         * 
         * \param script_path The path of the script to be loaded.
         */
    void LoadModule(const std::string& script_path);

    /**
         * \brief Calls the callable attribute "name" using the provided
         * arguments.
         * 
         * This function might throw a ERRORCODE if there is
         * an error when calling the function.
         * 
         * \param name The name of the attribute to be called.
         * \param args The arguments which will be used when calling the
         * attribute.
         * \return Python::Object containing the result of the function.
         */
    template <typename... Args>
    Object CallFunction(const std::string& name, const Args&... args)
    {
        //ensure global interpreter lock is aquired
        PyGILState_STATE state = PyGILState_Ensure();
        ON_SCOPE_EXIT([&] {PyGILState_Release(state); });
        Object func = GetAttr(name);
        ASSERT_ALLWAYS(PyCallable_Check(func.Get()), name << " is not a callable function!");
        // Create the tuple argument
        Object tup = PyTuple_New(sizeof...(args));
        add_tuple_vars(tup, args...);
        Object result = PyObject_CallObject(func.Get(), tup.Get());
        PropagatePyError();
        return result;
    }
    //
    //    /**
    //         * \brief Calls a callable attribute using no arguments.
    //         *
    //         * This function might throw a ERRORCODE if there is
    //         * an error when calling the function.
    //         *
    //         * \sa Python::Object::CallFunction.
    //         * \param name The name of the callable attribute to be executed.
    //         * \return Python::Object containing the result of the function.
    //         */
    Object CallFunction(const std::string& name);
    //
    //    /**
    //         * \brief Finds and returns the attribute named "name".
    //         *
    //         * This function might throw a ERRORCODE if an error
    //         * is encountered while fetching the attribute.
    //         *
    //         * \param name The name of the attribute to be returned.
    //         * \return Python::Object representing the attribute.
    //         */
    Object GetAttr(const std::string& name);
    //
    //    /**
    //         * \brief Checks whether this object contains a certain attribute.
    //         *
    //         * \param name The name of the attribute to be searched.
    //         * \return bool indicating whether the attribute is defined.
    //         */
    bool HasAttr(const std::string& name);

    /**
         * \brief Returns the internal PyObject*.
         * 
         * No reference increment is performed on the PyObject* before
         * returning it, so any DECREF applied to it without INCREF'ing
         * it will cause undefined behaviour.
         * \return The PyObject* which this Object is representing.
         */

private:
    // Variadic template method to add items to a tuple
    template <typename First, typename... Rest>
    void add_tuple_vars(Object& tup, const First& head, const Rest&... tail)
    {
        add_tuple_var(
            tup,
            PyTuple_Size(tup.Get()) - sizeof...(tail)-1,
            head);
        add_tuple_vars(tup, tail...);
    }

    void add_tuple_vars(Object& tup, Object arg)
    {
        add_tuple_var(tup, PyTuple_Size(tup.Get()) - 1, arg);
    }

    // Base case for add_tuple_vars
    template <typename Arg>
    void add_tuple_vars(Object& tup, const Arg& arg)
    {
        add_tuple_var(tup,
                      PyTuple_Size(tup.Get()) - 1, CastToPy(arg));
    }

    // Adds a Object to the tuple object
    void add_tuple_var(Object& tup, Py_ssize_t i, Object pobj)
    {
        PyTuple_SetItem(tup.Get(), i, pobj.Get(NewRef));
    }

    // Adds a PyObject* to the tuple object
    template <class T>
    void add_tuple_var(Object& tup, Py_ssize_t i,
                       const T& data)
    {
        Object item = CastToPy(data);
        PyTuple_SetItem(tup.Get(), i, item.Get(NewRef));
    }
};
};

#endif // PYWRAPPER_H
