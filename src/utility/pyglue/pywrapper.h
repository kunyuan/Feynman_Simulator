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
#include "object.h"

namespace Python {
class AnyObject : public Object {
public:
    AnyObject(Ownership ownership = Received);
    AnyObject(PyObject* ptr, Ownership ownership = Received);
    AnyObject(const AnyObject&, Ownership ownership = Received);
    template <typename T>
    AnyObject(T value)
    {
        _OwnerShip = Received;
        _PyPtr = CastToPyObject(value);
    }
    template <typename T>
    T As()
    {
        T value;
        if (Python::Convert(_PyPtr, value))
            ERRORCODEABORT(ERR_VALUE_INVALID, "Fail to convert PyObject!");
        return value;
    }

private:
    AnyObject& operator=(const AnyObject& a);
};

class ModuleObject : Object {
public:
    /**
         * \brief Constructs a default python object
         */
    ModuleObject();

    /**
         * \brief Constructs a python object from a PyObject pointer.
         * 
         * This Object takes(steals) ownership of the PyObject* argument. That 
         * means no Py_INCREF is performed on it. 
         * \param obj The pointer from which to construct this Object.
         */
    Object(PyObject* obj, Ownership ownership = TakeOver);
    template <typename T>
    Object(T obj, Ownership ownership = TakeOver)
        : py_obj(make_pyshared(CastToPyObject(obj), ownership))
    {
        _owenership = ownership;
    }
    ~Object()
    {
        //        LOG_INFO("Object ref: " << py_obj.use_count());
        //        LOG_INFO("pyObject ref: " << py_obj.get()->ob_refcnt);
    }

    void Print();
    std::string PrettyString();
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
        pyunique_ptr func(load_function(name));
        // Create the tuple argument
        pyunique_ptr tup(PyTuple_New(sizeof...(args)));
        add_tuple_vars(tup, args...);
        // Call our object
        PyObject* ret(PyObject_CallObject(func.get(), tup.get()));
        MakeSureNoError();
        return { ret };
    }

    /**
         * \brief Calls a callable attribute using no arguments.
         * 
         * This function might throw a ERRORCODE if there is
         * an error when calling the function.
         * 
         * \sa Python::Object::CallFunction.
         * \param name The name of the callable attribute to be executed.
         * \return Python::Object containing the result of the function.
         */
    Object CallFunction(const std::string& name);

    /**
         * \brief Finds and returns the attribute named "name".
         * 
         * This function might throw a ERRORCODE if an error
         * is encountered while fetching the attribute.
         * 
         * \param name The name of the attribute to be returned.
         * \return Python::Object representing the attribute.
         */
    Object GetAttr(const std::string& name);

    /**
         * \brief Checks whether this object contains a certain attribute.
         * 
         * \param name The name of the attribute to be searched.
         * \return bool indicating whether the attribute is defined.
         */
    bool HasAttr(const std::string& name);

    /**
         * \brief Returns the internal PyObject*.
         * 
         * No reference increment is performed on the PyObject* before
         * returning it, so any DECREF applied to it without INCREF'ing
         * it will cause undefined behaviour.
         * \return The PyObject* which this Object is representing.
         */
    PyObject* get() const { return py_obj.get(); }

    template <class T>
    bool Convert(T& param)
    {
        return Python::Convert(py_obj.get(), param);
    }
    template <class T>
    T As()
    {
        T param;
        if (Python::Convert(py_obj.get(), param))
            ERRORCODEABORT(ERR_VALUE_INVALID, "Fail to convert PyObject!");
        return param;
    }

    /**
         * \brief Constructs a Python::Object from a script file.
         * 
         * The returned Object will be the representation of the loaded
         * script. If any errors are encountered while loading this 
         * script, an ERRORCODE is thrown.
         * 
         * \param script_path The path of the script to be loaded.
         * \return Object representing the loaded script.
         */
    void LoadScript(const std::string& script_path);
    /**
         * \brief Constructs a Python::Object from a script string.
         * 
         * The returned Object will be the evaluation of the
         * script. If any errors are encountered while loading this 
         * script, an ERRORCODE is thrown.
         * 
         * \param script The string of the script to be evaluated.
         * \return Object representing the evaluated script.
         */
    void EvalScript(const std::string& script);
    void LoadModule(const std::string& module);

    size_t use_count() const { return py_obj.use_count(); }

private:
    typedef std::shared_ptr<PyObject> pyshared_ptr;
    pyshared_ptr py_obj;
    Ownership _owenership;
    void _deleter(PyObject* obj);

    PyObject* load_function(const std::string& name);

    pyshared_ptr make_pyshared(PyObject* obj, Ownership ownership);
    void reset_pyshared(PyObject* obj, Ownership ownership);

    // Variadic template method to add items to a tuple
    template <typename First, typename... Rest>
    void add_tuple_vars(pyunique_ptr& tup, const First& head, const Rest&... tail)
    {
        add_tuple_var(
            tup,
            PyTuple_Size(tup.get()) - sizeof...(tail)-1,
            head);
        add_tuple_vars(tup, tail...);
    }

    void add_tuple_vars(pyunique_ptr& tup, PyObject* arg)
    {
        add_tuple_var(tup, PyTuple_Size(tup.get()) - 1, arg);
    }

    // Base case for add_tuple_vars
    template <typename Arg>
    void add_tuple_vars(pyunique_ptr& tup, const Arg& arg)
    {
        add_tuple_var(tup,
                      PyTuple_Size(tup.get()) - 1, CastToPyObject(arg));
    }

    // Adds a PyObject* to the tuple object
    void add_tuple_var(pyunique_ptr& tup, Py_ssize_t i, PyObject* pobj)
    {
        PyTuple_SetItem(tup.get(), i, pobj);
    }

    // Adds a PyObject* to the tuple object
    template <class T>
    void add_tuple_var(pyunique_ptr& tup, Py_ssize_t i,
                       const T& data)
    {
        PyTuple_SetItem(tup.get(), i, CastToPyObject(data));
    }
};
};

#endif // PYWRAPPER_H
