
cp src/HODLR_ID_302.f temp/HODLR_ID_302.f
cp src/LUBackSolve.f temp/LUBackSolve.f
cp src/SMW_ZGESV2.f temp/SMW_ZGESV2.f

cp test/printTree.f temp/printTree.f
cp test/generateHODLR.f temp/generateHODLR.f
cp test/prini.f temp/prini.f
cp test/HODLR_test_302.f temp/HODLR_test_302.f
cp id_lib.a temp/id_lib.a


cd temp

gfortran -o HODLR_test HODLR_test_302.f HODLR_ID_302.f SMW_ZGESV2.f LUBackSolve.f generateHODLR.f printTree.f prini.f id_lib.a -llapack -lblas

./HODLR_test

rm HODLR_test HODLR_test_302.f HODLR_ID_302.f SMW_ZGESV2.f LUBackSolve.f readMatrix.f generateHODLR.f printTree.f prini.f id_lib.a
