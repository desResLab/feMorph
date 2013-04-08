#include "femModel.h"

femModel::femModel()
{
}

femModel::~femModel()
{
}

void femModel::ReadModelFromTFile(std::string fileName){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, is_any_of(" ,"), token_compress_on);
    std::string newstr = _copy("Hello World");
    
    if (boost::to_upper(tokenizedString[1])=='NODE'){
      // Read Nodes
    }else if (boost::to_upper(tokenizedString[1])=='TETRA10'){
      // Read Tetra 10 Element
    }else if (boost::to_upper(tokenizedString[1])=='PROP'){
      // Read Element Property
    }
  }
  
  // Close File
  infile.close();
}

