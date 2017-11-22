clc;
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O is_border_valsIMPORT.cpp
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O random_initIMPORT.cpp
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O populate_indices.cpp
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O SP_prop_init.cpp

cd optical_flow_celiu/mex/
mex -O Coarse2FineTwoFrames.cpp GaussianPyramid.cpp OpticalFlow.cpp
cd ../../

mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O localonly_move.cpp IMG.cpp NormalD.cpp SP.cpp;
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O local_move.cpp IMG.cpp NormalD.cpp SP.cpp;
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O switch_move.cpp IMG.cpp NormalD.cpp SP.cpp;
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O merge_move.cpp IMG.cpp NormalD.cpp SP.cpp;
mex -Igsl/include -Lgsl/lib -lgsl -lcblas -O split_move.cpp IMG.cpp NormalD.cpp SP.cpp;
