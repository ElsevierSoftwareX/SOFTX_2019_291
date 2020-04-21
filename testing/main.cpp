#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#define TST_NO_EXEC_OUT
#include "testing.hpp"

void test_convert(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    #### Convert test suite ####\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test conv_bin("meshtool convert -imsh=meshes/block/block.elem -omsh=/tmp/test_block.belem"),
                      conv_vtk("meshtool convert -imsh=meshes/block/block.elem -omsh=/tmp/test_block.vtk"),
                      conv_vtu("meshtool convert -imsh=meshes/block/block.elem -omsh=/tmp/test_block.vtu");

  conv_bin.add_check(new tst::bin_file_equal("equal binary elem test",
                     "convert/ref_block.belem", "/tmp/test_block.belem"));
  conv_bin.add_check(new tst::bin_file_equal("equal binary lon test",
                     "convert/ref_block.blon", "/tmp/test_block.blon"));
  conv_bin.add_check(new tst::bin_file_equal("equal binary pts test",
                     "convert/ref_block.bpts", "/tmp/test_block.bpts"));

  conv_vtk.add_check(new tst::bin_file_equal("equal binary file test",
                     "convert/ref_block.vtk", "/tmp/test_block.vtk"));

  conv_vtu.add_check(new tst::bin_file_equal("equal binary file test",
                     "convert/ref_block.vtu", "/tmp/test_block.vtu"));

  tst::execution_test rem_tmp("rm /tmp/test_block.b* /tmp/test_block.vtk /tmp/test_block.vtu");

  error += conv_bin.run();
  error += conv_vtk.run();
  error += conv_vtu.run();
  error += rem_tmp.run();

  if(error)
    sprintf(summary + strlen(summary), "Suite \"conversion tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"conversion tests\" succeeded with no errors.\n\n");
}


void test_extract(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    #### Extract test suite ####\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test surf("meshtool extract surface -msh=meshes/block/block.elem -surf=/tmp/block"),
                      submsh("meshtool extract mesh -msh=meshes/block/block.elem -tags=2 -submsh=/tmp/test_tag2.vtk"),
                      data("meshtool extract data -msh=meshes/block/block -submsh=extract/tag2 "
                           "-msh_data=extract/block.act.dat -submsh_data=/tmp/tag2.act.dat");

  surf.add_check(new tst::text_file_equal("equal .surf.vtx file",
                     "extract/block.surf.vtx", "/tmp/block.surf.vtx"));

  submsh.add_check(new tst::bin_file_equal("equal binary file",
                     "extract/tag2.vtk", "/tmp/test_tag2.vtk"));

  data.add_check(new tst::float_text_file_equal("equal .dat file",
                 "extract/tag2.act.dat", "/tmp/tag2.act.dat", 1e-6));

  tst::execution_test rem_tmp("rm /tmp/block.* /tmp/test_tag2.vtk /tmp/tag2.act.dat");

  error += int(surf.run());
  error += int(submsh.run());
  error += int(data.run());
  error += int(rem_tmp.run());

  if(error)
    sprintf(summary + strlen(summary), "Suite \"extraction tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"extraction tests\" succeeded with no errors.\n\n");
}

void test_insert(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    #### Insert test suite ####\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test submsh("meshtool insert submesh -msh=meshes/block/block.elem"
                             " -submsh=extract/tag2 -outmsh=/tmp/insert.block"),
                      data("meshtool insert data -msh=meshes/block/block -submsh=extract/tag2 "
                           "-msh_data=extract/block.act.dat -submsh_data=extract/tag2.act.dat -odat=/tmp/insert.dat");

  submsh.add_check(new tst::text_file_equal("equal .elem file",
                     "meshes/block/block.elem", "/tmp/insert.block.elem"));
  submsh.add_check(new tst::float_text_file_equal("equal .pts file",
                     "meshes/block/block.pts", "/tmp/insert.block.pts", 1e-6));
  submsh.add_check(new tst::float_text_file_equal("equal .lon file",
                     "meshes/block/block.lon", "/tmp/insert.block.lon", 1e-6));

  data.add_check(new tst::float_text_file_equal("equal .dat file",
                 "extract/block.act.dat", "/tmp/insert.dat", 1e-6));

  tst::execution_test rem_tmp("rm /tmp/insert.block.* /tmp/insert.dat");

  error += int(submsh.run());
  error += int(data.run());
  error += int(rem_tmp.run());

  if(error)
    sprintf(summary + strlen(summary), "Suite \"insertion tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"insertion tests\" succeeded with no errors.\n\n");
}

void test_resample(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    ###  Resample test suite ###\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test mesh_ref("meshtool resample mesh -msh=meshes/block/block.elem -outmsh=/tmp/block.ref.vtk -max=200 -postsmth=0"),
                      mesh_crs("OMP_NUM_THREADS=1 meshtool resample mesh -msh=resample/block.ref.vtk -outmsh=/tmp/block.crs.vtk -min=700 -surf_corr=0.95 -postsmth=0");

  mesh_ref.add_check(new tst::bin_file_equal("equal .vtk file",
                     "resample/block.ref.vtk", "/tmp/block.ref.vtk"));

  mesh_crs.add_check(new tst::bin_file_equal("equal .vtk file",
                     "resample/block.crs.vtk", "/tmp/block.crs.vtk"));

  tst::execution_test rem_tmp("rm /tmp/block.*");

  error += int(mesh_ref.run());
  error += int(mesh_crs.run());
  error += int(rem_tmp.run());

  if(error)
    sprintf(summary + strlen(summary), "Suite \"resample tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"resample tests\" succeeded with no errors.\n\n");
}

void test_smooth(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    ### Smoothing test suite ###\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test mesh_smth("OMP_NUM_THREADS=1 meshtool smooth mesh -msh=resample/block.ref.vtk -outmsh=/tmp/block.ref.smth.elem -smth=0.2 -iter=200 -tags=+ -thr=0.95 -edge=30");

  mesh_smth.add_check(new tst::text_file_equal("equal .elem file",
                     "smooth/block.ref.smth.elem", "/tmp/block.ref.smth.elem"));
  mesh_smth.add_check(new tst::float_text_file_equal("equal .pts file",
                     "smooth/block.ref.smth.pts", "/tmp/block.ref.smth.pts", 1e-3));
  mesh_smth.add_check(new tst::float_text_file_equal("equal .lon file",
                     "smooth/block.ref.smth.lon", "/tmp/block.ref.smth.lon", 1e-3));

  tst::execution_test rem_tmp("rm /tmp/block.ref.smth.*");

  error += int(mesh_smth.run());
  error += int(rem_tmp.run());

  if(error)
    sprintf(summary + strlen(summary), "Suite \"smoothing tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"smoothing tests\" succeeded with no errors.\n\n");
}


void test_indexsearch(char summary[]) {

  printf("\n\n    ############################\n");
  printf(    "    ### Index search test suite ###\n");
  printf(    "    ############################\n");

  int error = 0;

  tst::execution_test query_idxlist("meshtool query idxlist -msh=meshes/block/block -coord=indexsearch/searchfile.txt");

  query_idxlist.add_check(new tst::text_file_equal("equal output textfile",
                          "indexsearch/searchfile.txt.out.txt", "indexsearch/reference.results.txt"));

  tst::execution_test rem_tmp("rm indexsearch/searchfile.txt.out.txt");

  error += int(query_idxlist.run());
  error += int(rem_tmp.run());

  if(error)
    sprintf(summary + strlen(summary), "Suite \"indexsearch tests\" failed with %d errors.\n\n", error);
  else
    sprintf(summary + strlen(summary), "Suite \"indexsearch tests\" succeeded with no errors.\n\n");
}



int main(int argc, char** argv)
{
  const int bufflength = 100000;
  char summary[bufflength];
  sprintf(summary, "\n\n");

  test_convert(summary);
  test_extract(summary);
  test_insert(summary);
  test_resample(summary);
  test_smooth(summary);
  test_indexsearch(summary);

  printf("%s", summary);
  return 0;
}
