#include "Settings.hpp"
#include "Coord.hpp"
#include "Grid.hpp"
#include "Error.hpp"

using namespace std;

const char* Error::what() {
	return "Unknown error.";
}

const char* IndexError::what() {
	return "Index out of boundaries error.";
}

const char* NoDataError::what() {
	return "No required data available.";
}

const char* NoFieldError::what() {
	return "No required up-to-date field available.";
}

const char* ConfigError::what() {
	return "Error parsing configuration file.";
}

const char* NoSolutionError::what() {
	return "No numeric solution found.";
}
