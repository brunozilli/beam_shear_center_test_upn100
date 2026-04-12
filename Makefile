# Makefile for beam_shear_centre
# Fortran 77 fixed format with LAPACK
# MIT License
# Copyright (c) 2024 Bruno Zilli & DeepSeek

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

# Test executables
TESTS = $(TESTDIR)/test_shear_center

all: $(OBJECTS) $(TESTS)

# Rule for test_shear_center
$(TESTDIR)/test_shear_center: $(TESTDIR)/test_shear_center.o $(OBJECTS)
	$(FC) -o $@ $^ $(LDFLAGS)

# Generic compilation rule
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(SRCDIR)/*.o $(TESTDIR)/*.o $(TESTS)

# Phony targets
.PHONY: all clean test test-all

# Run test on HEB200
test: $(TESTDIR)/test_shear_center
	./$(TESTDIR)/test_shear_center meshes/HEB200_mm.unv

# Run test on all meshes
test-all: $(TESTDIR)/test_shear_center
	@echo "========================================="
	@echo "  BEAM SHEAR CENTRE - MULTI-MESH TEST"
	@echo "========================================="
	@echo ""
	@for mesh in meshes/*.unv; do \
		echo "-----------------------------------------"; \
		echo "TESTING: $$(basename $$mesh)"; \
		echo "-----------------------------------------"; \
		./$(TESTDIR)/test_shear_center $$mesh; \
		echo ""; \
	done
	@echo "========================================="
	@echo "  TEST COMPLETED"
	@echo "========================================="
