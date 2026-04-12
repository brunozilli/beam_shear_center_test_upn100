#!/bin/bash

echo "========================================="
echo "  BEAM SHEAR CENTRE - MULTI-MESH TEST"
echo "========================================="
echo ""

for mesh in meshes/*.unv; do
    echo "-----------------------------------------"
    echo "TESTING: $(basename "$mesh")"
    echo "-----------------------------------------"
    ./test/test_shear_center "$mesh"
    echo ""
done

echo "========================================="
echo "  TEST COMPLETED"
echo "========================================="
