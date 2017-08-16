#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include "defaults.hpp"
#include "filterstream.hpp"

namespace unipeak {

  using namespace std;
  using namespace boost;
  using namespace boost::iostreams;

  void InStream::close () {
    if (input_) {
      input_->reset();
      delete input_;
      input_ = 0;
    }
  }

  void InStream::open (const string& fnameArg) {
    close();
    input_ = new filtering_istream;
    if (fnameArg == STDIN_STRING) {
      input_->push(cin);
      fname_ = "standard input stream";
    } else {
      if (ends_with(fnameArg, GZIP_SUFFIX)) {
        input_->push(gzip_decompressor());
      } else if (ends_with(fnameArg, BZIP2_SUFFIX)) {
        input_->push(bzip2_decompressor());
      }
      file_source in (fnameArg.c_str());
      if (! in.is_open()) {
        cerr << "error: could not read " << fnameArg << endl << endl;
        exit(1);
      }
      input_->push(in);
      fname_ = fnameArg;
    }
  }

  InStream::InStream () : input_(0), lineNo_(0) {}

  InStream::InStream (const string& fnameArg) : input_(0), lineNo_(0) {
    open(fnameArg);
  }

  InStream::~InStream () {
    close();
  }

  bool InStream::good () const {
    if (! input_) return false; else return(input_->good());
  }

  string InStream::readLine () {
    assert(good());
    string output;
    getline(*input_, output);
    ++lineNo_;
    return output;
  }

  const Counter& InStream::getLineNo () const {
    return lineNo_;
  }


  void OutStream::close () {
    if (output_) {
      output_->reset();
      delete output_;
      output_ = 0;
    }
  }

  void OutStream::open (const string& fnameArg) {
    close();
    output_ = new filtering_ostream;
    if (fnameArg == STDOUT_STRING) {
      output_->push(cout);
      fname_ = "standard output stream";
    } else {
      if (ends_with(fnameArg, GZIP_SUFFIX)) {
        output_->push(gzip_compressor());
      } else if (ends_with(fnameArg, BZIP2_SUFFIX)) {
        output_->push(bzip2_compressor());
      }
      file_sink out (fnameArg.c_str());
      if (! out.is_open()) {
        cerr << "error: could not write " << fnameArg << endl << endl;
        exit(1);
      }
      output_->push(out);
      fname_ = fnameArg;
    }
  }

  OutStream::OutStream () : output_(0) {}

  OutStream::OutStream (const string& fname) : output_(0) {
    open(fname);
  }

  OutStream::~OutStream () {
    close();
  }

  bool OutStream::good () const {
    return(output_->good());
  }

  void OutStream::write (const string& input) {
    assert(good());
    (*output_) << input;
  }


  OutStream& operator<< (OutStream& output, const string& input) {
    output.write(input);
    return output;
  }

  OutStream& operator<< (OutStream& output, const double& input) {
    output.write(lexical_cast<string>(input));
    return output;
  }

}
