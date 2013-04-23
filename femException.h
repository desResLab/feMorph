#ifndef FEMEXCEPTION_H
#define FEMEXCEPTION_H

#include <string>
#include <exception>

class femException: public std::exception{
public:
  // Constructor and Destructor
  femException(const char* m):msg(m){};
  virtual ~femException() throw(){};
    // Member Functions
      virtual const char* what() const throw() {return msg.c_str();}
protected:
  std::string msg;
};

#endif // FEMEXCEPTION_H
