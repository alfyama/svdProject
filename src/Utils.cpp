#include "Utils.h"

void readMatrixCsv(std::string fileName, std::vector<float>& mat, int &m, int &n){
   std::ifstream file(fileName);
   std::string line;

   int rows = 0;
   int cols = 0;

   while (std::getline(file, line)) {
       std::istringstream iss(line);
       float num;

       while (iss >> num) {
           mat.push_back(num);
           cols++;
       }
       rows++;
       if(rows == 1) cols /= rows;
   }
   m = rows;
   n = cols;
}


void readMatrixCsv(std::string fileName, std::vector<double>& mat, int &m, int &n){
   std::ifstream file(fileName);
   std::string line;

   int rows = 0;
   int cols = 0;

   while (std::getline(file, line)) {
       std::istringstream iss(line);
       double num;

       while (iss >> num) {
           mat.push_back(num);
           cols++;
       }
       rows++;
       if(rows == 1) cols /= rows;
   }
   m = rows;
   n = cols;
}


void readMatrixCsv(std::string fileName, std::vector<long double>& mat, int &m, int &n){
   std::ifstream file(fileName);
   std::string line;

   int rows = 0;
   int cols = 0;

   while (std::getline(file, line)) {
       std::istringstream iss(line);
       long double num;

       while (iss >> num) {
           mat.push_back(num);
           cols++;
       }
       rows++;
       if(rows == 1) cols /= rows;
   }
   m = rows;
   n = cols;
}
