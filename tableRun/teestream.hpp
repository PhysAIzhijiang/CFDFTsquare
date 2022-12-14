// This file is under public domain.
#ifndef TEESTREAM_HPP
#define TEESTREAM_HPP
#include <streambuf>
#include <ostream>

/*
 * The code is copied from
 * http://wordaligned.org/articles/cpp-streambufs
 *
 * Usage example:

#include <fstream>
#include <iostream>
#include <teestream>

int main()
{
    std::ofstream log("hello-world.log");
    teestream tee(std::cout, log);
    tee << "Hello, world!\n";
    return 0;
}
 */

template <typename char_type,
          typename traits = std::char_traits<char_type> >
class basic_teebuf:
    public std::basic_streambuf<char_type, traits>
{
public:
    typedef typename traits::int_type int_type;
    
    basic_teebuf(std::basic_streambuf<char_type, traits> * sb1,
                 std::basic_streambuf<char_type, traits> * sb2)
      : sb1(sb1)
      , sb2(sb2)
    {
    }
    
private:    
    virtual int sync()
    {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
    
    virtual int_type overflow(int_type c)
    {
        int_type const eof = traits::eof();
        
        if (traits::eq_int_type(c, eof))
        {
            return traits::not_eof(c);
        }
        else
        {
            char_type const ch = traits::to_char_type(c);
            int_type const r1 = sb1->sputc(ch);
            int_type const r2 = sb2->sputc(ch);
            
            return
                traits::eq_int_type(r1, eof) ||
                traits::eq_int_type(r2, eof) ? eof : c;
        }
    }
    
private:
    std::basic_streambuf<char_type, traits> * sb1;
    std::basic_streambuf<char_type, traits> * sb2;
};


typedef basic_teebuf<char> teebuf;

class teestream : public std::ostream
{
public:
    // Construct an ostream which tees output to the supplied
    // ostreams.
    teestream(std::ostream & o1, std::ostream & o2);
private:
    teebuf tbuf;
};

teestream::teestream(std::ostream & o1, std::ostream & o2)
  : std::ostream(&tbuf)
  , tbuf(o1.rdbuf(), o2.rdbuf())
{
}

#endif
