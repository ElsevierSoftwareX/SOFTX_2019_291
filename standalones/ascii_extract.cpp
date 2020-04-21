/**
* @file ascii_extract.cpp
* @brief Standalone binary for reading components from ascii vector data files.
* @author Aurel Neic
* @version 
* @date 2017-09-12
*/
#include<stdio.h>

#include<iostream>
#include<string>

#include "mt_utils.h"


template<class T>
void parse_nodes(const std::string & nodes, mt_vector<T> & nod_vec)
{
  std::set<T> nod_set;
  mt_vector<std::string> oper_strings;
  split_string(nodes, ',' , oper_strings);

  for(size_t op = 0; op < oper_strings.size(); op++)
  {
    mt_vector<std::string> range_strings;
    split_string(oper_strings[op], '-', range_strings);

    if(range_strings.size() == 2)
    {
      T start = atoi(range_strings[0].c_str());
      T end   = atoi(range_strings[1].c_str());
      for(T i = start; i<end; i++)
        nod_set.insert(i);
    }
    else
      nod_set.insert(atoi(range_strings[0].c_str()));
  }

  nod_vec.resize(nod_set.size());

  size_t idx = 0;
  for(auto it = nod_set.begin(); it != nod_set.end(); ++it)
    nod_vec[idx++] = *it;
}



int main(int argc, char** argv)
{
  if(argc < 5)
  {
    fprintf(stderr, "Error: Wrong usage!\n");
    fprintf(stderr, "Use:\n");
    fprintf(stderr, "%s nodefile <nodes file> <data in> <data out>\n", argv[0]);
    fprintf(stderr, "Or:\n");
    fprintf(stderr, "%s nodes <nodes> <data in> <data out>\n", argv[0]);
    exit(1);
  }

  std::string mode        = argv[1];
  std::string input_file  = argv[3];
  std::string output_file = argv[4];
  mt_vector<mt_int> nodes;

  if(mode.compare("nodefile") == 0)
  {
    std::string nodes_file = argv[2];
    read_vector_ascii(nodes, nodes_file, false);
  }
  else if(mode.compare("nodes") == 0)
  {
    std::string nodes_string = argv[2];
    parse_nodes(nodes_string, nodes);
  }
  else {
    std::cerr << "Error paring extracting mode." << std::endl;
  }

  std::cout << "Reading input data .." << std::endl;
  mt_vector<mt_real> input;
  read_vector_ascii(input, input_file, true);

  size_t num_nodes_in  = input.size();
  size_t num_nodes_out = nodes.size();
  std::cout << "Input data size: " << num_nodes_in << " output data size: " << num_nodes_out << std::endl;

  std::cout << "Mapping data .." << std::endl;
  mt_vector<mt_real> output(num_nodes_out);

  for(size_t i=0; i < num_nodes_out; i++) {
    output[i] = input[nodes[i]];
  }

  std::cout << "Writing output data .." << std::endl;
  write_vector_ascii(output, output_file);
  return 0;
}




