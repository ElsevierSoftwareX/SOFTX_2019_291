include ../my_switches.def
include ../switches.def

CXXOPTS += -I../utils -I../modes
LIB = ../utils/libmtutils.a ../modes/libmtmodes.a

ifdef MT_SILENT
CXX := @$(CXX)
endif

default: correspondance ascii_extract samesideofplane retag simplify_tag_topology

dynpts_tools: convert_dynpts combine_dynpts transfer refine_time_igb fillup_dynpts deactivate_igb_tags

all: enclose retag simplify_tag_topology correspondance ascii_extract vtxToMesh fibmap convert_dynpts combine_dynpts transfer refine_time_igb fillup_dynpts samesideofplane convert_igb_to_flowvc volumefractions convert_obj_to_igb deactivate_igb_tags

enclose: enclose.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o enclose enclose.cpp $(LIB)

retag: retag.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o retag retag.cpp $(LIB)

simplify_tag_topology: simplify_tag_topology.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o simplify_tag_topology simplify_tag_topology.cpp $(LIB)

tara_msh_retag: tara_msh_retag.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o tara_msh_retag tara_msh_retag.cpp $(LIB)

reaction_force: reaction_force.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o reaction_force reaction_force.cpp $(LIB)

correspondance: correspondance.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o correspondance correspondance.cpp $(LIB)

ascii_extract: ascii_extract.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o ascii_extract ascii_extract.cpp $(LIB)

vtxToMesh: vtxToMesh.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o vtxToMesh vtxToMesh.cpp $(LIB)

fibmap: fibmap.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o fibmap fibmap.cpp $(LIB)

spctm_data: spctm_data.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o spctm_data spctm_data.cpp $(LIB)

convert_dynpts: convert_xdynpts_to_displdynpts.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o convert_dynpts convert_xdynpts_to_displdynpts.cpp $(LIB)

combine_dynpts: combine_dynpts_in_basemesh.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o combine_dynpts combine_dynpts_in_basemesh.cpp $(LIB)

partition: partition.cpp ../utils/kdpart.hpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o partition partition.cpp $(LIB)

transfer: transfer_data_elasticity_fluid.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o transfer transfer_data_elasticity_fluid.cpp $(LIB)

refine_time_igb: refine_time_igb.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o refine_time_igb refine_time_igb.cpp $(LIB)

smooth_mesh_transient: smooth_mesh_transient.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o smooth_mesh_transient smooth_mesh_transient.cpp $(LIB)

fillup_dynpts: fillup_dynpts.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o fillup_dynpts fillup_dynpts.cpp $(LIB)

change_dynpts_header: change_dynpts_header.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o change_dynpts_header change_dynpts_header.cpp $(LIB)

calculate_bad_elems_transient: calculate_bad_elems_transient.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o calculate_bad_elems_transient calculate_bad_elems_transient.cpp $(LIB)

samesideofplane: samesideofplane.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o samesideofplane samesideofplane.cpp $(LIB)

convert_igb_to_flowvc: convert_igb_to_flowvc.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o convert_igb_to_flowvc convert_igb_to_flowvc.cpp $(LIB)

volumefractions: volumefractions.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o volumefractions volumefractions.cpp $(LIB)

convert_obj_to_igb: convert_obj_to_igb.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o convert_obj_to_igb convert_obj_to_igb.cpp $(LIB)

deactivate_igb_tags: deactivate_igb_based_on_tags.cpp
	@echo "* compiling:" $<
	$(CXX) $(CXXOPTS) -o deactivate_igb_tags deactivate_igb_based_on_tags.cpp $(LIB)


clean:
	rm -f enclose correspondance fibmap ascii_extract retag simplify_tag_topology \
tara_msh_retag vtxToMesh convert_dynpts combine_dynpts transfer refine_time_igb \
smooth_mesh_transient fillup_dynpts change_dynpts_header reaction_force samesideofplane \
convert_igb_to_flowvc spctm_data partition volumefractions convert_obj_to_igb deactivate_igb_tags

