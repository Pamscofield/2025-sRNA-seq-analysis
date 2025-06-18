
python ~/Downloads/software/SparK-master/SparK.py \
-gs yes \
-cf CEA17_merge_minus.bdg CEA17_dox_merge_minus.bdg CEA17_pabA_IRT_merge_minus.bdg \
-tf CEA17_merge_plus.bdg CEA17_dox_merge_plus.bdg CEA17_pabA_IRT_merge_plus.bdg \
-tg 1 2 3 \
-cg 1 2 3 \
-gl CEA17 CEA17_dox CEA17_pabA_IRT \
-l - + \
-pr chrDS499602:177871-178419 \
-o pabA \
-f 0072B2 E69F00 \
-pt sine

python ~/Downloads/software/SparK-master/SparK.py \
-gs yes \
-cf CEA17_merge_minus.bdg CEA17_dox_merge_minus.bdg dclA_pksp_IRT_merge_minus.bdg dclB_pksp_IRT_merge_minus.bdg \
-tf CEA17_merge_plus.bdg CEA17_dox_merge_plus.bdg dclA_pksp_IRT_merge_plus.bdg dclB_pksp_IRT_merge_plus.bdg \
-tg 1 2 3 4 \
-cg 1 2 3 4 \
-gl CEA17 CEA17_dox dclA_pksp_IRT dclB_pksp_IRT \
-l - + \
-pr chrDS499595:4646340-4646839 \
-o pksp \
-f 0072B2 E69F00 \
-pt sine
