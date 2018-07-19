#ifndef _ERROR_HPP
#define _ERROR_HPP
#include <iostream>

class Error: public std::exception {
	public:
	const char* what();
};

class IndexError: public Error {
	public:
	const char* what();
};

class ConfigError: public Error {
	public:
	const char *what();
};
	
class NoDataError: public Error {
	public:
	const char *what();
};

class NoFieldError: public NoDataError {
	public:
	const char *what();
};

class NoSolutionError: public Error {
	public:
	const char *what();
};



#endif
