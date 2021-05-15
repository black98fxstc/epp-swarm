//
// Created by Jonathan Ebrahimian on 10/16/20.
//

#ifndef UTILS_FILELOADER_H
#define UTILS_FILELOADER_H
#include <string>
#include <fstream>
#include <vector>

class FileLoader {
private:
    //selected info
    bool headers;
    bool display;
    int rows;
    int cols;
    char delimiter;

    int dimensions;
    int dataType;
    std::string fileName;
    std::ifstream inFile;
    int storageType;
    bool verbose;

    //data
    std::vector<std::vector<int>> int2DV;
    std::vector<std::vector<std::string>> string2DV;
    std::vector<std::vector<double>> double2DV;
    std::vector<int> int1DV;
    std::vector<double> double1DV;
    std::vector<std::string> string1DV;

    int ** int2DA;
    std::string ** string2DA;
    double ** double2DA;
    int * int1DA;
    double * double1DA;
    std::string * string1DA;



    //Funcs
    void skipHeader();
    void parse();
    void displayData();
    void allocateMem();
    void addData(std::vector<std::string> &,int);
    void deAllocateMem();
    void breakLine(std::string,std::vector<std::string>&,bool);
    void print2DArray(int);
    void print1DArray(int);
    void print2DVector(int);
    void print1DVector(int);





public:
    FileLoader();
    enum fileTypes{
        CSV = 0,
        TXT
    };
    enum storageTypes{
        VECTOR = 0,
        ARRAY
    };
//    enum fileTypes{
//
//    };
    enum type {
        INT = 0,
        DOUBLE,
        STRING
    };
    enum delimiters{
        SPACE = 0,
        COMMA,
        SEMICOLON,
        COLON,
        BANG,
        PIPE
    };
    void selectDimensions(int);
    void selectDataType(int);
    void selectDataType(std::string);
    void selectStorageType(int);
    void selectStorageType(std::string);
    void selectVerbose(bool);
    void selectCols(int);
    void selectRows(int);
    void selectDelimiter(char);
    void selectDelimiter(int);
    void selectDisplay(bool);
    void selectHeaders(bool);
    void execute();
    void selectDisplayCheck(bool);
    void selectFileName(std::string);
    void checkRun();
    void checkRun(int);
    void readME();
    void refresh();



    void * get();








};


#endif //UTILS_FILELOADER_H
