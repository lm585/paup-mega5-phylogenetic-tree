/*
$ head  cucumber-snp-7chrs.miss0d05-minAllele2.25Ksnp.transpose | cut -f 1-10
accession       Chr1_9765       Chr1_11951      Chr1_11989      Chr1_27336      Chr1_33527      Chr1_41394      Chr1_42563      Chr1_58201      Chr1_72147
CG1697_Ea       C       C       C       G       G       A       A       G       G
CG5801_Ea       C       C       C       G       G       A       A       G       G
CG1811_Ea       C       C       C       G       K       A       A       G       G
CG3085_Ea       C       C       C       G       G       A       A       G       G
CG3087_Ea       C       C       C       A       G       A       A       G       G
CG4081_Ea       C       C       C       G       G       A       A       G       G
CG4182_Ea       C       C       C       A       T       A       A       G       G
CG4210_Ea       C       C       C       A       G       A       A       G       G
CG1778_Ea       C       C       C       A       T       A       A       G       G

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
 if(argc != 4)
 {
  cerr << argv[0] << "  snp-genotype-table-transposed    1  output-file" << endl;
  cerr << "snp-table:   only capital A/C/G/T   K/M/R/S/W/Y" << endl;
  cerr << "snp-table:  tab delimited table. See the head of the C++ program" << endl;
  cerr << "1: ignore the 1st line of the snp-table" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, string> base2int;
 int ignoreLine;

 base2int["A"] = "11";
 base2int["C"] = "22";
 base2int["G"] = "33";
 base2int["T"] = "44";
 base2int["K"] = "34"; // K,k ==> {GT}
 base2int["M"] = "12"; //M,m ==> {AC}
 base2int["R"] = "13"; //R,r ==> {AG}
 base2int["S"] = "23"; //S,s ==> {CG}
 base2int["W"] = "14"; //W,w ==> {AT}
 base2int["Y"] = "24"; //Y,y ==> {CT}
 

 ignoreLine = atoi(argv[2]);
 readFile(argv[1], input);
 writeFile(argv[3], output);

 for(int i = 0; i < ignoreLine; i++)
 {
  getline(input, line);
 }
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  for(int i = 0; i < 2; i++)
  {
   output << lineFields[0] << '\t'; 
   for(int j = 1; j < lineFields.size(); j++)
   {
    string s;
    s = lineFields[j];
    if(s == "-")
    {
     output << -9;
    }
    else if(base2int.count(s) == 0)
    {
     output << -7;
    }
    else  //ACGTKMR...
    {
     output << base2int[s][i];
    }
    if(j == lineFields.size() - 1)
      output << endl;
    else
      output << '\t';
   }
  }
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


