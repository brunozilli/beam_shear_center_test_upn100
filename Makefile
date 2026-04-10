# Makefile for beam_shear_centre
# Fortran 77 fixed format with LAPACK

FC = gfortran
FFLAGS = -ffixed-form -Wall -O2 -std=legacy
LDFLAGS = -llapack -lblas

SRCDIR = src
TESTDIR = test

SOURCES = $(SRCDIR)/read_section_mesh_unv.f \
          $(SRCDIR)/mesh_checker.f \
          $(SRCDIR)/compute_section_properties.f \
          $(SRCDIR)/shear_center.f \
          $(SRCDIR)/build_D_matrix.f \
          $(SRCDIR)/section_database.f

OBJECTS = $(SOURCES:.f=.o)
TESTS = $(TESTDIR)/test_shear_center

all: $(OBJECTS) $(TESTS)

$(TESTDIR)/test_shear_center: $(TESTDIR)/test_shear_center.o $(OBJECTS)
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(TESTDIR)/*.o $(TESTS)

test_shear_center: $(TESTDIR)/test_shear_center

.PHONY: all clean test_shear_center
