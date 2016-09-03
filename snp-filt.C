/*
input:
Chr1_770        C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C - C C C C C C C C C C C C C C C Y C C C C C C C C C C C C C C C C C C - C C C C T T C C C C C C C C C C - C C C - C C C C C C C C C C C C C C C C C C C C C C C C C C C C
Chr1_780        C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C Y Y C C T C C C C C C C C C C C C C C - C C C - C C C C C C C C C C C C C C C C C C C C C C C C C C C C
Chr1_1155       T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T - T T T Y Y Y T T T T T T T T T T T T T T T T T T T T T T T T T Y Y T T T T T T T T T T T T T T T T T T - T T T - T T T T T T T T T T T T T T T T T T T T T T T T T T T T
Chr1_13791      C C C C T C C C C C C C Y C C C C C C Y C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C Y C C C C C C C C C C C C C C C Y C C C C C C C Y C T C C T T T C C C Y C C C T C Y C Y C C T T T T T T T T T T T T T T T T T T T
Chr1_14376      C C Y C T Y T T T Y C C T T Y T C C C Y C C C C Y C T Y T C - C T C C C C C Y C C C C C C C C C C C Y C C C C C C C C C C C C - C C Y C C C C C C T Y T T T Y T T T T C C T T C C T C T C Y T C T T T T T T T T T T T T T T T T T T T
Chr1_14407      C C Y C T Y T T T Y C C T T Y T C C C Y C C C C Y C T Y T C C C T C C C C C Y C C C C C C C Y C C C Y C C C C C C C C C C C C C C C Y C C C T T T T Y T T T Y T T T T C C T T C C T C - C C T T T T T T T T T T T T T T T T T T T T T

input:  max missing rate (0.05)
input: minor allele detected in at least 2 accessions

op: output file

 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 5)
 {
  cerr << argv[0] << " snp_table  0.05    2    kept-snp-file-name" << endl;
  cerr << "snp_table: delimited by space or tab" << endl;
  cerr << "max missing rate (0.05)" << endl;
  cerr << "minor allele detected in at least 2 accessions" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 double missRate;
 int  minorAllele;

 missRate = atof(argv[2]);
 minorAllele = atoi(argv[3]);
 writeFile(argv[4], output);
 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, "\t\n ", line);
   int na, nc, ng, nt, miss, numAllele;
   na = 0;
   nc = 0;
   ng = 0;
   nt = 0;
   miss = 0;
   numAllele = 0;
  for(int i = 1; i < lineFields.size(); i++) //skip first column
  {
   if(lineFields[i] == "A")
     na++;
   else if(lineFields[i] == "C")
     nc++;
   else if(lineFields[i] == "G")
     ng++;
   else if(lineFields[i] == "T")
     nt++;
   else if(lineFields[i] == "-")
     miss++;
   else;
  }
  if(na >= minorAllele) numAllele++;
  if(nc >= minorAllele) numAllele++;
  if(ng >= minorAllele) numAllele++;
  if(nt >= minorAllele) numAllele++;
  if(numAllele >= 2 && miss * 1.0 / (lineFields.size() - 1.0) < missRate)
    output << line << endl;
   
  getline(input, line);
 }
 input.close();
 
 output.close();
 return 0;
}

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out);
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

void getFieldContent(vector<string> & fields, char del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del == line[i])
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}

void getFieldContent2(vector<string> & fields, string del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del.find(line[i]) != string::npos)
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   if(len > 0)
   {
    str = line.substr(pos[i] + 1, len);
    fields.push_back(str);
   }
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}


