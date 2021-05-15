//
// Created by Jonathan Ebrahimian on 10/16/20.
//

#include "FileLoader.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <limits>

FileLoader::FileLoader() {
    rows = std::numeric_limits<int>::max();
    cols = std::numeric_limits<int>::max();
    headers = false;
    display = false;
    delimiter = ',';
    dimensions = 2;
    dataType = type::INT;
    storageType = storageTypes::VECTOR;
    fileName = "";

    int2DA = nullptr;
    string2DA = nullptr;
    double2DA = nullptr;
    int1DA = nullptr;
    double1DA = nullptr;
    string1DA = nullptr;

}

void FileLoader::selectDimensions(int dimensionsIn) {
    dimensions = dimensionsIn;
}

void FileLoader::selectDataType(int typeIn) {
    if(typeIn == type::INT){
        dataType = type::INT;
    }else if(typeIn == type::DOUBLE){
        dataType = type::DOUBLE;
    }else if(typeIn == type::STRING){
        dataType = type::STRING;
    }else {
        std::cout << "Invalid cast type selected" << std::endl;
    }
}

void FileLoader::selectDataType(std::string typeIn) {
    std::transform(typeIn.begin(), typeIn.end(), typeIn.begin(), ::tolower);
    if(typeIn == "int"){
        dataType = type::INT;
    }else if(typeIn == "double"){
        dataType = type::DOUBLE;
    }else if(typeIn == "string"){
        dataType = type::STRING;
    }else {
        std::cout << "Invalid cast type selected" << std::endl;
    }
}

void FileLoader::selectStorageType(int storageInput){

    if(storageInput == storageTypes::VECTOR){
        storageType = storageTypes::VECTOR;
    }else if(storageInput == storageTypes::ARRAY){
        storageType = storageTypes::ARRAY;
        if(verbose){
            std::cout << "Remember to set row and cols." << std::endl;
        }
    }else {
        std::cout << "Invalid cast type selected" << std::endl;
    }

}
void FileLoader::selectStorageType(std::string inputStorageType){
    std::transform(inputStorageType.begin(), inputStorageType.end(), inputStorageType.begin(), ::tolower);
    if(inputStorageType == "vector"){
        storageType = storageTypes::VECTOR;
    }else if(inputStorageType == "array"){
        storageType = storageTypes::ARRAY;
        if(verbose){
            std::cout << "remember to set row and cols." << std::endl;
        }
    }else {
        std::cout << "Invalid cast type selected" << std::endl;
    }
}

void FileLoader::selectVerbose(bool verboseIn) {
    verbose = verboseIn;
}

void FileLoader::selectCols(int colsIn) {
    cols = colsIn;
}

void FileLoader::selectRows(int rowsIn) {
    rows = rowsIn;
}

void FileLoader::selectDelimiter(char delimiterIn) {
    delimiter = delimiterIn;
}

void FileLoader::selectDelimiter(int delimiterIn) {
    if(delimiterIn == delimiters::SPACE){
        delimiter = ' ';
    }else if(delimiterIn == delimiters::COMMA){
        delimiter = ',';
    }else if(delimiterIn == delimiters::PIPE){
        delimiter = '|';
    }else if(delimiterIn == delimiters::COLON){
        delimiter = ':';
    }else if(delimiterIn == delimiters::BANG){
        delimiter = '!';
    }else if(delimiterIn == delimiters::SEMICOLON){
        delimiter = ';';
    }else{
        std::cout << "Invalid delimiter selected" << std::endl;
    }
}

void FileLoader::selectDisplay(bool) {

}

void FileLoader::selectHeaders(bool headersIn) {
    headers = headersIn;
}

void FileLoader::execute() {
    if(cols == std::numeric_limits<int>::max() && storageType == storageTypes::ARRAY && dimensions == 2){
        std::cout << "Must set custom col number for array storage type" << std::endl;
        exit(0);
    }else if(rows == std::numeric_limits<int>::max() && storageType == storageTypes::ARRAY){
        std::cout << "Must set custom row number for array storage type"  << std::endl;
        exit(0);
    }
    inFile.open(fileName);
    if(!inFile.is_open()){
        std::cout << "Could not open file:  " << fileName << std::endl;
        std::cout << "Make sure you have given the correct path to the file." << std::endl;
        exit(0);
    }

    deAllocateMem();
    allocateMem();

    if(headers){
        skipHeader();
    }
    parse();
    if(display){
        displayData();
    }

    inFile.close();
}

void FileLoader::selectDisplayCheck(bool displayIn) {
    display = displayIn;
}

void FileLoader::selectFileName(std::string fileNameIn) {
    fileName = fileNameIn;
    if(storageType == storageTypes::VECTOR){
        rows = std::numeric_limits<int>::max();
        cols = std::numeric_limits<int>::max();
    }

}

void FileLoader::skipHeader() {
    std::string tempLine;
    getline(inFile,tempLine);
}

void FileLoader::parse() {
    int numRowsParse = 0;
    std::string tempLine;
    std::vector<std::string> stringVec;
    bool flag = true;
    if(cols != std::numeric_limits<int>::max() && dimensions == 2){
        stringVec.resize(cols);
        flag = false;
    }


    while(!inFile.eof() && numRowsParse < rows){
        getline(inFile,tempLine);
        breakLine(tempLine,stringVec,flag);
        if(dimensions == 2) {
            if (storageType == storageTypes::VECTOR) {
                if(dataType == type::INT){
                    int2DV.push_back(std::vector<int>());
                }else if(dataType == type::DOUBLE){
                    double2DV.push_back(std::vector<double>());
                }else if(dataType == type::STRING){
                    string2DV.push_back(stringVec);
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }
        }
        addData(stringVec,numRowsParse);
        numRowsParse += 1;
    }

    //setting row and col value from dynamic read
    if(flag){
        if(dimensions == 2) {
            if (storageType == storageTypes::VECTOR) {
                if(dataType == type::INT){
                    rows = int2DV.size();
                    //cols = int2DV[0].size();
                }else if(dataType == type::DOUBLE){
                    rows = double2DV.size();
                    //cols = double2DV[0].size();
                }else if(dataType == type::STRING){
                    rows = string2DV.size();
                    //cols = string2DV[0].size();
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }
        }else if(dimensions == 1){
            if (storageType == storageTypes::VECTOR) {
                if(dataType == type::INT){
                    rows = int1DV.size();
                }else if(dataType == type::DOUBLE){
                    rows = double1DV.size();
                }else if(dataType == type::STRING){
                    rows = string1DV.size();
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }
        }else{
            std::cout << "Invalid dimmension selected." << std::endl;
        }
    }
}
void FileLoader::addData(std::vector<std::string> & vec,const int numRowsParse) {
    for(int x = 0; x < vec.size();x++){
        if(dimensions == 1){
            if(storageType == storageTypes::VECTOR){
                if(dataType == type::INT){
                    int1DV.push_back(std::stoi(vec[x]));
                }else if(dataType == type::DOUBLE){
                    double1DV.push_back(std::stod(vec[x]));
                }else if(dataType == type::STRING){
                    string1DV.push_back(vec[x]);
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }else if(storageType == storageTypes::ARRAY){
                if(dataType == type::INT){
                    int1DA[numRowsParse] = std::stoi(vec[x]);
                }else if(dataType == type::DOUBLE){
                    double1DA[numRowsParse] = std::stod(vec[x]);
                }else if(dataType == type::STRING){
                    string1DA[numRowsParse] = vec[x];
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }else{
                std::cout << "Invalid storage type selected." << std::endl;
            }
        }else if(dimensions == 2){
            //String 2D vec not added because it is taken care of above
            if(storageType == storageTypes::VECTOR){
                if(dataType == type::INT){
                    int2DV[numRowsParse].push_back(std::stoi(vec[x]));
                }else if(dataType == type::DOUBLE){
                    double2DV[numRowsParse].push_back(std::stod(vec[x]));
                }else if(dataType == type::STRING){

                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }else if(storageType == storageTypes::ARRAY){
                if(dataType == type::INT){
                    int2DA[numRowsParse][x] = std::stoi(vec[x]);
                }else if(dataType == type::DOUBLE){
                    double2DA[numRowsParse][x] = std::stod(vec[x]);
                }else if(dataType == type::STRING){
                    string2DA[numRowsParse][x] = vec[x];
                }else{
                    std::cout << "Invalid data type selected." << std::endl;
                }
            }else{
                std::cout << "Invalid storage type selected." << std::endl;
            }
        }else{
            std::cout << dimensions << "  is an invalid dimension" << std::endl;
            exit(0);
        }

    }
}
void FileLoader::breakLine(std::string line, std::vector<std::string> & vec,bool reset){
    if(reset){
        vec.clear();
    }
    int initialSize = vec.size();
    int colCount = 0;
    std::stringstream ss(line);
    std::string segment;
    while(std::getline(ss, segment, delimiter))
    {
        if(initialSize == 0){
            vec.push_back(segment);
        }else{
            if(colCount < cols){
                vec[colCount] = segment;
            }else{
                break;
            }

        }
        colCount += 1;
    }
    if(vec[vec.size()-1] == ""){
        vec.pop_back();
    }
}

void FileLoader::displayData() {

}


void FileLoader::allocateMem() {
    if(dimensions == 1){
        if(storageType == storageTypes::VECTOR){
//            if(dataType == type::INT){
//                int1DV.resize(rows);
//            }else if(dataType == type::DOUBLE){
//                double1DV.resize(rows);
//            }else if(dataType == type::STRING){
//                string1DV.resize(rows);
//            }else{
//                std::cout << "Invalid data type selected." << std::endl;
//            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                int1DA = new int[rows];
            }else if(dataType == type::DOUBLE){
                double1DA = new double[rows];
            }else if(dataType == type::STRING){
                string1DA = new std::string[rows];
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else if(dimensions == 2){
        if(storageType == storageTypes::VECTOR){
//            if(dataType == type::INT){
//                int2DV = std::vector<std::vector<int>>(rows, std::vector<int>(cols, 0));
//            }else if(dataType == type::DOUBLE){
//                double2DV = std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0));
//            }else if(dataType == type::STRING){
//                string2DV = std::vector<std::vector<std::string>>(rows, std::vector<std::string>(cols, 0));
//            }else{
//                std::cout << "Invalid data type selected." << std::endl;
//            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                int2DA = new int*[rows];
                for(int i = 0; i < rows; ++i)
                    int2DA[i] = new int[cols];
            }else if(dataType == type::DOUBLE){
                double2DA = new double*[rows];
                for(int i = 0; i < rows; ++i)
                    double2DA[i] = new double[cols];
            }else if(dataType == type::STRING){
                string2DA = new std::string*[rows];
                for(int i = 0; i < rows; ++i)
                    string2DA[i] = new std::string[cols];
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }

    }else{
        std::cout << dimensions << "  is an invalid dimension" << std::endl;
        exit(0);
    }
}

void FileLoader::deAllocateMem() {
    if(verbose && storageType == ARRAY){
        std::cout << "Make sure to get your data before reading in a new file to avoid a guarenteed memory leak." << std::endl;
    }
    if(dimensions == 1){
        if(storageType == storageTypes::VECTOR){
            if(dataType == type::INT){
                int1DV.clear();
            }else if(dataType == type::DOUBLE){
                double1DV.clear();
            }else if(dataType == type::STRING){
                string1DV.clear();
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                int1DA = nullptr;
            }else if(dataType == type::DOUBLE){
                double1DA = nullptr;
            }else if(dataType == type::STRING){
                string1DA = nullptr;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else if(dimensions == 2){
        if(storageType == storageTypes::VECTOR){
            if(dataType == type::INT){
                int2DV.clear();
            }else if(dataType == type::DOUBLE){
                double2DV.clear();
            }else if(dataType == type::STRING){
                string2DV.clear();
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                int2DA = nullptr;
            }else if(dataType == type::DOUBLE){
                double2DA = nullptr;
            }else if(dataType == type::STRING){
                string2DA = nullptr;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }

    }else{
        std::cout << dimensions << "  is an invalid dimension" << std::endl;
        exit(0);
    }
}

void *FileLoader::get() {
    if(dimensions == 1){
        if(storageType == storageTypes::VECTOR){
            if(dataType == type::INT){
                return &int1DV;
            }else if(dataType == type::DOUBLE){
                return &double1DV;
            }else if(dataType == type::STRING){
                return &string1DV;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                return int1DA;
            }else if(dataType == type::DOUBLE){
                return double1DA;
            }else if(dataType == type::STRING){
                return string1DA;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else if(dimensions == 2){
        if(storageType == storageTypes::VECTOR){
            if(dataType == type::INT){
                return &int2DV;
            }else if(dataType == type::DOUBLE){
                return &double2DV;
            }else if(dataType == type::STRING){
                return &string2DV;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else if(storageType == storageTypes::ARRAY){
            if(dataType == type::INT){
                return int2DA;
            }else if(dataType == type::DOUBLE){
                return double2DA;
            }else if(dataType == type::STRING){
                return string2DA;
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else{
        std::cout << "Select a dimension before getting data" << std::endl;
        exit(0);
    }
}

void FileLoader::readME() {
    std::cout << "VARS and their defaults: " << std::endl;
    std::cout << "rows = INF" << std::endl;
    std::cout << "cols = INF" << std::endl;
    std::cout << "headers = false" << std::endl;
    std::cout << "display = false" << std::endl;
    std::cout << "dimensions = 2" << std::endl;
    std::cout << "dataType = INT" << std::endl;
    std::cout << "storageType = VECTOR" << std::endl;
    std::cout << "filename = \"\"" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Presets will read a 2D file of ints and store it in an array." << std::endl;
    std::cout << "Remember to type cast the void pointer returned when getting data." << std::endl;
}

void FileLoader::print2DArray(int numRows) {
    if(numRows == -1){
        numRows = rows;
    }
    for(int x = 0; x < rows; x++){
        for(int y = 0; y < cols;y++){
            if(dataType == type::INT){
                std::cout << int2DA[x][y] << " ";
            }else if(dataType == type::DOUBLE){
                std::cout << double2DA[x][y] << " ";
            }else if(dataType == type::STRING){
                std::cout << string2DA[x][y] << " ";
            }else{
                std::cout << "Invalid data type selected." << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

void FileLoader::print1DArray(int numRows) {
    if(numRows == -1){
        numRows = rows;
    }
    for(int y = 0; y < numRows ;y++){
        if(dataType == type::INT){
            std::cout << int1DA[y] << std::endl;
        }else if(dataType == type::DOUBLE){
            std::cout << double1DA[y] <<std::endl;
        }else if(dataType == type::STRING){
            std::cout << string1DA[y] << std::endl;
        }else{
            std::cout << "Invalid data type selected." << std::endl;
        }
    }
}

void FileLoader::print2DVector(int numRows) {
    if(numRows == -1){
        numRows = rows;
    }

    if(dataType == type::INT){
        for(int x = 0; x < rows; x++){
            for(int y = 0; y < int2DV[x].size();y++){
                std::cout << int2DV[x][y] << " ";
            }
            std::cout << std::endl;
        }
    }else if(dataType == type::DOUBLE){
        for(int x = 0; x < rows; x++){
            for(int y = 0; y < double2DV[x].size();y++){
                std::cout << double2DV[x][y] << " ";
            }
            std::cout << std::endl;
        }
    }else if(dataType == type::STRING){
        for(int x = 0; x < rows; x++){
            for(int y = 0; y < string2DV[x].size();y++){
                std::cout << string2DV[x][y] << " ";
            }
            std::cout << std::endl;
        }
    }else{
        std::cout << "Invalid data type selected." << std::endl;
    }

}

void FileLoader::print1DVector(int numRows) {
    if(numRows == -1){
        numRows = rows;
    }
    for(int y = 0; y < numRows ;y++){
        if(dataType == type::INT){
            std::cout << int1DV[y] << std::endl;
        }else if(dataType == type::DOUBLE){
            std::cout << double1DV[y] <<std::endl;
        }else if(dataType == type::STRING){
            std::cout << string1DV[y] << std::endl;
        }else{
            std::cout << "Invalid data type selected." << std::endl;
        }
    }
}

void FileLoader::checkRun() {
    if(dimensions == 1){
        if(storageType == storageTypes::VECTOR){
            print1DVector(-1);
        }else if(storageType == storageTypes::ARRAY){
            print1DArray(-1);
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else if(dimensions == 2){
        if(storageType == storageTypes::VECTOR){
            print2DVector(-1);
        }else if(storageType == storageTypes::ARRAY){
            print2DArray(-1);
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else{
        std::cout << "Select a dimension before getting data" << std::endl;
        exit(0);
    }
}

void FileLoader::checkRun(int numRows) {
    if(dimensions == 1){
        if(storageType == storageTypes::VECTOR){
            print1DVector(numRows);
        }else if(storageType == storageTypes::ARRAY){
            print1DArray(numRows);
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else if(dimensions == 2){
        if(storageType == storageTypes::VECTOR){
            print2DVector(numRows);
        }else if(storageType == storageTypes::ARRAY){
            print2DArray(numRows);
        }else{
            std::cout << "Invalid storage type selected." << std::endl;
        }
    }else{
        std::cout << "Select a dimension before getting data" << std::endl;
        exit(0);
    }
}

void FileLoader::refresh(){
    rows = std::numeric_limits<int>::max();
    cols = std::numeric_limits<int>::max();
    headers = false;
    display = false;
    delimiter = ',';
    dimensions = 2;
    dataType = type::INT;
    storageType = storageTypes::VECTOR;
    fileName = "";
}
