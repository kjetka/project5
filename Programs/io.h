#ifndef IO_H
#define IO_H
#include <fstream>
class System;
using std::ofstream;

class IO{
private:
    ofstream file;
public:
    IO(const char *filename);
    ~IO();          //spr ~???

    void saveState(System &system);
    void open(const char *filename);
    void close();

};
#endif
