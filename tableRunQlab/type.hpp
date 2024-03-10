// This file is under public domain.
#ifndef TYPE_HPP
#define TYPE_HPP
/*
 * https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
 */
#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}


/*
int& foo_lref();
int&& foo_rref();
int foo_value();

int
main()
{
    int i = 0;
    const int ci = 0;
    std::cout << "decltype(i) is " << type_name<decltype(i)>() << '\n';
    std::cout << "decltype((i)) is " << type_name<decltype((i))>() << '\n';
    std::cout << "decltype(ci) is " << type_name<decltype(ci)>() << '\n';
    std::cout << "decltype((ci)) is " << type_name<decltype((ci))>() << '\n';
    std::cout << "decltype(static_cast<int&>(i)) is " << type_name<decltype(static_cast<int&>(i))>() << '\n';
    std::cout << "decltype(static_cast<int&&>(i)) is " << type_name<decltype(static_cast<int&&>(i))>() << '\n';
    std::cout << "decltype(static_cast<int>(i)) is " << type_name<decltype(static_cast<int>(i))>() << '\n';
    std::cout << "decltype(foo_lref()) is " << type_name<decltype(foo_lref())>() << '\n';
    std::cout << "decltype(foo_rref()) is " << type_name<decltype(foo_rref())>() << '\n';
    std::cout << "decltype(foo_value()) is " << type_name<decltype(foo_value())>() << '\n';
}
 */
#endif
