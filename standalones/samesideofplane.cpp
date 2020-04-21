/**
* @file samesideofplane.cpp
* @brief Standalone binary. A plane is define through 3 nodes. Based on a pts-file and a
*                           reference node the code will decide on if the test pt is located
*                           on the same side as the reference point or on the opposite.
* @author Anton Prassl
* @version
* @date 2018-04-10
*/
#include<stdio.h>
#include<iostream>
#include<string>

#include "mt_utils.h"
#include "mt_modes_base.h"

#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

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


int sameSideOfPlane(const mt_point<mt_real> pa,
                    const mt_point<mt_real> pb,
                    const mt_point<mt_real> pc,
                    const mt_point<mt_real> pref,
                    const mt_point<mt_real> ptest)
{
    mt_real value;

    mt_point<mt_real>        p01 = pb - pa;
    mt_point<mt_real>        p02 = pc - pa;
    mt_point<mt_real> inw_normal = p01.crossProd(p02);

    value       = inw_normal.scaProd(pref-pa);
    inw_normal *= SIGN(value);

    if (inw_normal.scaProd(ptest-pa) < 0.)
      return 0; // point on different side
    else
      return 1;  // point on same side
}


int main(int argc, char** argv)
{
  struct timeval t1, t2;

  if(argc < 16)
  {
    fprintf(stderr, "Error: Wrong usage!\n");
    fprintf(stderr, "Use:\n");
    fprintf(stderr, "%s ptsfile datfile const ptA ptB ptC pRef \n", argv[0]);
    fprintf(stderr, "Or:\n");
    fprintf(stderr, "%s nodes <nodes> <data in> <data out>\n", argv[0]);
    exit(1);
  }
  // start global timer
  gettimeofday(&t1, NULL);

  std::string ptsFile = argv[1];
  std::string datFile = argv[2];
  int value           = atoi(argv[3]);
  mt_point<mt_real> ptA( atof(argv[4]), atof(argv[5]), atof(argv[6]) );
  mt_point<mt_real> ptB( atof(argv[7]), atof(argv[8]), atof(argv[9]) );
  mt_point<mt_real> ptC( atof(argv[10]), atof(argv[11]), atof(argv[12]) );
  mt_point<mt_real> ptRef( atof(argv[13]), atof(argv[14]), atof(argv[15]) );

  size_t numpts = readNumPoints(ptsFile);
  mt_vector<mt_real> nodes(3*numpts);
  mt_vector<mt_int> data(numpts,0);

  std::cout << "Reading " << ptsFile << std::endl;
  readPoints(nodes, ptsFile);
  //binary_write(nodes.begin(), nodes.end(), "myptsFile.pts");

  if (file_exists(datFile)) {
    std::cout << "Reading " << datFile << std::endl;
    read_vector_ascii(data, datFile);
  } else
    std::cout << "Creating " << datFile << std::endl;

  PROGRESS<size_t> progress(numpts, "Evaluating test nodes: ");

  for(size_t i=0; i<numpts; i++)
  {
    mt_point<mt_real> ptTest(nodes.data()+3*i);

    progress.next();
    if (sameSideOfPlane(ptA, ptB, ptC, ptRef, ptTest))
       data[i] += (mt_int)value;
  }
  progress.finish();


  std::cout << "Writing output data .." << std::endl;
  write_vtx(data, datFile, false);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}




