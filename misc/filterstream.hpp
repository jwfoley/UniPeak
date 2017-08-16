#ifndef FILTERSTREAM_H
#define FILTERSTREAM_H

#include <cstdlib>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>

#include "defaults.hpp"

namespace unipeak {

  using namespace std;
  using namespace boost::iostreams;

  class InStream {
    public:
      InStream ();
      InStream (const string&);
      ~InStream ();
      void open (const string&);
      void close ();
      bool good () const;
      string readLine ();
      const Counter& getLineNo () const;
      
    protected:
      filtering_istream* input_;
      string fname_;
      Counter lineNo_;
  };

  class OutStream {
    public:
      OutStream ();
      OutStream (const string&);
      ~OutStream ();
      void open (const string&);
      void close ();
      bool good () const;
      void write (const string&);
    
    protected:
      filtering_ostream* output_;
      string fname_;
  };

  OutStream& operator<< (OutStream&, const string&);
  OutStream& operator<< (OutStream&, const double&);

}

#endif
