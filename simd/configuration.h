#ifndef ___CONFIGURATION___
#define ___CONFIGURATION___

#include <fstream>

class Configuration {
public:
	static unsigned int cv_k_fold;
	static unsigned int classNumber;
	static std::ofstream communication;
};

#endif
