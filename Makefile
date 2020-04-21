include my_switches.def
include switches.def

OBJ_UTILS = \
utils/data_structs.o \
utils/geometry_utils.o \
utils/io_utils.o \
utils/string_utils.o \
utils/mesh_quality.o \
utils/mesh_utils.o \
utils/mesh_refining.o \
utils/mesh_smoothing.o \
utils/mesh_generation.o \
utils/topology_utils.o \
utils/mmg_utils.o \
utils/stellar_utils.o \
utils/vtk_utils.o \
utils/vcflow_utils.o \
utils/netgen_utils.o \
utils/ensight_utils.o \
utils/mt_igb.o \
utils/kdtree.o \
utils/tetgen/tetgen.o \
utils/tetgen/predicates.o

OBJ_MODES = \
modes/clean_mode.o \
modes/convert_mode.o \
modes/collect_mode.o \
modes/extract_mode.o \
modes/generate_mode.o \
modes/insert_mode.o \
modes/interpolate_mode.o \
modes/itk_mode.o \
modes/map_mode.o \
modes/merge_mode.o \
modes/mt_modes_base.o \
modes/query_mode.o \
modes/resample_mode.o \
modes/reindex_mode.o \
modes/restore_mode.o \
modes/smooth_mode.o \
modes/split_mode.o \
modes/transform_mode.o

DEPS = $(OBJ_UTILS:%.o=%.d) $(OBJ_MODES:%.o=%.d)

LIB_UTILS = utils/libmtutils.a
LIB_MODES = modes/libmtmodes.a
LIBS = -Lutils -Lmodes -lmtmodes -lmtutils

ifdef MT_ADDONS
CXXOPTS += -DMT_ADDONS -Imeshtool_addons/src
LIBS += -Lmeshtool_addons -l addons
endif

default: parbuild

all: parbuild stdl

parbuild:
	make $(MKFLG) meshtool_deps
	make $(MKSLFLG) meshtool

stdl:
	cd standalones; make $(MKFLG) all

doc: README.html Advanced.Workflows.html

my_switches.def: switches.def
	@echo "* generating "$@
	@echo "#MT_DEBUG = 1         # Compile in debug mode" > my_switches.def
	@echo "#MT_STATIC = 1        # Try to link statically" >> my_switches.def
	@echo "MT_OPENMP = 1        # Use OpenMP" >> my_switches.def
	@echo "#MT_ADDONS = 1        # Include addons" >> my_switches.def
	@echo "MT_SILENT = 1        # Dont show build command" >> my_switches.def
	@echo "MT_LARGE_TYPES = 1   # Whether to use the (long,double) or (int,float) for (mt_int,mt_real)" >> my_switches.def
	@echo "MT_CC_ENV = gnu      # choose between intel, gnu, clang" >> my_switches.def
	@echo "PAR_MAKE = 8         # choose the number of processes to use for parallel make" >> my_switches.def
	@echo "#MT_SILENT_PRG = 1    # use silent progress instead of a bar" >> my_switches.def
	@echo "#WINDOWS_BUILD = 1    # turn on when building on windows" >> my_switches.def

utils/%.o: utils/%.cpp
	@echo "* compiling:" $<
	$(CXX) -c -o $@ $(CXXOPTS) $<

utils/tetgen/%.o: utils/tetgen/%.cxx
	@echo "* compiling:" $<
	$(CXX) -c -o $@ $(CXXOPTS) $<

modes/%.o: modes/%.cpp
	@echo "* compiling:" $<
	$(CXX) -c -o $@ -Iutils $(CXXOPTS) $<

%.html : %.md
	python md_to_html.py $< -o $@

$(LIB_UTILS): $(OBJ_UTILS)
	@echo "* creating archive:" $@
	ar cr $(LIB_UTILS) utils/*.o utils/tetgen/*.o
	ranlib $(LIB_UTILS)

$(LIB_MODES): $(OBJ_MODES)
	@echo "* creating archive:" $@
	ar cr $(LIB_MODES) modes/*.o
	ranlib $(LIB_MODES)

meshtool_deps: $(LIB_MODES) $(LIB_UTILS)

meshtool: main.cpp meshtool_deps
ifdef MT_ADDONS
	@echo "* building addons .."
	@(cd meshtool_addons && make $(MKFLG) libaddons.a)
endif
	@echo "* linking:" $@
	$(CXX) $(CXXOPTS) -Imodes -Iutils -o meshtool main.cpp $(LIBS)

doxygen:
	doxygen doxydoc/mt_doc.doxygen

clean:
	rm -f meshtool
	rm -f utils/*.[od] utils/tetgen/*.[od] utils/*.a
	rm -f modes/*.[od] modesmodes/*.a
	rm -rf doxydoc/html
	cd standalones && make clean
ifdef MT_ADDONS
	cd meshtool_addons && make clean
endif

status:
	git status
ifdef MT_ADDONS
	cd meshtool_addons && git status
endif

update:
	git pull
ifdef MT_ADDONS
	cd meshtool_addons && git pull
endif

-include $(DEPS)
