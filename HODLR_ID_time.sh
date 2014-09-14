
cp src/HODLR_ID_301.f temp/HODLR_ID_301.f

cp test/printTree.f temp/printTree.f
cp test/generateHODLR.f temp/generateHODLR.f
cp test/prini.f temp/prini.f
cp test/HODLR_timeTest.f temp/HODLR_timeTest.f
cp id_lib.a temp/id_lib.a


cd temp

gfortran -o HODLR_test HODLR_timeTest.f HODLR_ID_301.f generateHODLR.f printTree.f prini.f id_lib.a -llapack -lblas
./HODLR_test

rm HODLR_test HODLR_timeTest.f HODLR_ID_301.f readMatrix.f generateHODLR.f printTree.f prini.f id_lib.a
