
cp src/HODLR_ID_302.f temp/HODLR_ID_302.f
cp src/HODLR_ID_det.f temp/HODLR_ID_det.f
cp src/LUBackSolve.f temp/LUBackSolve.f
cp src/SMW_ZGESV2.f temp/SMW_ZGESV2.f
cp src/DGETRF_2.f temp/DGETRF_2.f

cp test/printTree.f temp/printTree.f
cp test/generateHODLR.f temp/generateHODLR.f
cp test/prini.f temp/prini.f
cp test/HODLR_test_det.f temp/HODLR_test_det.f
cp id_lib.a temp/id_lib.a


cd temp

gfortran -o HODLR_test_det HODLR_test_det.f HODLR_ID_302.f HODLR_ID_det.f DGETRF_2.f SMW_ZGESV2.f LUBackSolve.f generateHODLR.f printTree.f prini.f id_lib.a -llapack -lblas

./HODLR_test_det

rm HODLR_test_det HODLR_test_det.f HODLR_ID_302.f HODLR_ID_det.f DGETRF_2.f SMW_ZGESV2.f LUBackSolve.f generateHODLR.f printTree.f prini.f id_lib.a
