CC = gfortran
TESTS_FOLDER = tests/
TEST = bspline_action
TEST_SRC = $(TEST).f90
EXE = $(TEST)
MODS_SRC = bspline_sub_module.f90 bspline_oo_module.f90 bspline_module.f90 bspline.f90
MODS = $(MODS_SRC:.f90=.mod)
all: $(EXE)
	

$(EXE): $(TESTS_FOLDER)$(TEST_SRC) $(MODS) $(MODS_SRC)
	$(CC) $(CFLAGS) $(TESTS_FOLDER)$(TEST_SRC) $(MODS_SRC) -o $(EXE)

%.mod: %.f90
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.mod *.o $(EXES)
