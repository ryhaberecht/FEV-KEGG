"""
Question
--------
Which reaction edges (substrate -> reaction ID -> product) are present in eco01100, but not in a combination of all non-overview metabolic pathways?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Calculate difference (eco01100 - allNonOverviewPathways). Do not subtract nodes, to keep all surviving edges.
- Print reaction edges (substrate -> reaction ID -> product).

Result
------

::

    877 results
    C00002 -> R02014 -> C00131
    C00003 -> R00189 -> C00857
    C00003 -> R00257 -> C00857
    C00006 -> R00104 -> C00003
    C00010 -> R00130 -> C00882
    C00011 -> R00214 -> C00149
    C00011 -> R00216 -> C00149
    C00011 -> R00345 -> C00036
    C00014 -> R00149 -> C00169
    C00014 -> R01221 -> C00037
    C00016 -> R00160 -> C00061
    C00016 -> R00161 -> C00061
    C00018 -> R00173 -> C00250
    C00018 -> R00174 -> C00250
    C00018 -> R00277 -> C00647
    C00018 -> R00278 -> C00627
    C00019 -> R00177 -> C00073
    C00019 -> R10404 -> C00021
    C00020 -> R00185 -> C00212
    C00020 -> R01083 -> C03794
    C00021 -> R04858 -> C00019
    C00021 -> R10404 -> C00019
    C00022 -> R00196 -> C00186
    C00022 -> R00200 -> C00074
    C00022 -> R00206 -> C00074
    C00022 -> R00209 -> C00024
    C00022 -> R00220 -> C00065
    C00022 -> R00369 -> C00041
    C00022 -> R00396 -> C00041
    C00022 -> R00400 -> C00041
    C00022 -> R00703 -> C00186
    C00022 -> R01447 -> C00186
    C00022 -> R02320 -> C00074
    C00022 -> R03105 -> C00957
    C00024 -> R00209 -> C00022
    C00024 -> R00229 -> C00033
    C00024 -> R00352 -> C00158
    C00025 -> R00093 -> C00026
    C00025 -> R00114 -> C00026
    C00025 -> R00243 -> C00026
    C00025 -> R00355 -> C00049
    C00025 -> R00411 -> C05931
    C00026 -> R00243 -> C00025
    C00026 -> R00267 -> C00311
    C00026 -> R00709 -> C00311
    C00029 -> R00287 -> C00103
    C00031 -> R00010 -> C01083
    C00031 -> R00028 -> C00208
    C00031 -> R00305 -> C00198
    C00031 -> R00306 -> C00185
    C00031 -> R00801 -> C00089
    C00031 -> R01555 -> C00208
    C00031 -> R02727 -> C01083
    C00031 -> R05196 -> C00208
    C00031 -> R07147 -> C00198
    C00032 -> R00310 -> C02191
    C00033 -> R00229 -> C00024
    C00033 -> R00320 -> C00227
    C00033 -> R00710 -> C00084
    C00036 -> R00345 -> C00011
    C00036 -> R00345 -> C00074
    C00036 -> R00346 -> C00074
    C00036 -> R00352 -> C00158
    C00036 -> R00357 -> C00049
    C00036 -> R00361 -> C00149
    C00036 -> R00431 -> C00074
    C00036 -> R00726 -> C00074
    C00037 -> R00899 -> C01419
    C00039 -> R00377 -> C00458
    C00039 -> R00378 -> C00459
    C00041 -> R00369 -> C00022
    C00041 -> R00396 -> C00022
    C00041 -> R00400 -> C00022
    C00041 -> R03599 -> C05688
    C00042 -> R00432 -> C00091
    C00042 -> R00713 -> C00232
    C00042 -> R00714 -> C00232
    C00042 -> R00727 -> C00091
    C00042 -> R10343 -> C00091
    C00043 -> R00414 -> C00645
    C00043 -> R00416 -> C04501
    C00044 -> R02020 -> C00286
    C00046 -> R00442 -> C00063
    C00046 -> R00443 -> C00075
    C00047 -> R00451 -> C00680
    C00048 -> R00475 -> C00160
    C00048 -> R00479 -> C00311
    C00048 -> R00717 -> C00160
    C00049 -> R00355 -> C00025
    C00049 -> R00485 -> C00152
    C00049 -> R07407 -> C05840
    C00049 -> R07410 -> C05840
    C00051 -> R00497 -> C00037
    C00051 -> R00497 -> C00669
    C00051 -> R01918 -> C05730
    C00052 -> R00502 -> C00446
    C00053 -> R00509 -> C00224
    C00055 -> R00513 -> C00475
    C00055 -> R00516 -> C00475
    C00055 -> R00517 -> C00475
    C00055 -> R00962 -> C00475
    C00055 -> R01548 -> C00475
    C00055 -> R02091 -> C00475
    C00055 -> R02096 -> C00475
    C00055 -> R02371 -> C00475
    C00055 -> R02372 -> C00475
    C00058 -> R00519 -> C00080
    C00059 -> R00530 -> C00224
    C00059 -> R00531 -> C00224
    C00061 -> R00160 -> C00016
    C00061 -> R00549 -> C00255
    C00061 -> R00550 -> C00255
    C00061 -> R05705 -> C01847
    C00061 -> R08574 -> C00255
    C00062 -> R01086 -> C03406
    C00063 -> R00571 -> C00075
    C00063 -> R00573 -> C00075
    C00063 -> R02022 -> C00458
    C00064 -> R00253 -> C00014
    C00064 -> R00253 -> C00025
    C00065 -> R00582 -> C01005
    C00068 -> R00616 -> C03028
    C00068 -> R00617 -> C01081
    C00068 -> R00618 -> C03028
    C00073 -> R00648 -> C01180
    C00073 -> R00650 -> C00155
    C00073 -> R00946 -> C00155
    C00073 -> R02821 -> C00155
    C00073 -> R04405 -> C00155
    C00073 -> R07396 -> C01180
    C00074 -> R00199 -> C00022
    C00074 -> R00206 -> C00022
    C00074 -> R00346 -> C00036
    C00074 -> R00431 -> C00036
    C00074 -> R00726 -> C00036
    C00077 -> R00669 -> C00437
    C00077 -> R02282 -> C00437
    C00078 -> R00674 -> C00463
    C00078 -> R02722 -> C00065
    C00079 -> R00688 -> C00166
    C00079 -> R00689 -> C00166
    C00079 -> R00692 -> C00166
    C00080 -> R00519 -> C00058
    C00082 -> R00729 -> C01179
    C00082 -> R09254 -> C01179
    C00084 -> R00710 -> C00033
    C00084 -> R00746 -> C00469
    C00084 -> R05198 -> C00469
    C00084 -> R09127 -> C00469
    C00085 -> R00762 -> C00354
    C00085 -> R01067 -> C00279
    C00091 -> R00432 -> C00042
    C00091 -> R00727 -> C00042
    C00091 -> R10343 -> C00042
    C00093 -> R00841 -> C00116
    C00093 -> R00847 -> C00116
    C00094 -> R00859 -> C00283
    C00094 -> R00861 -> C00283
    C00094 -> R02021 -> C00053
    C00095 -> R00801 -> C00089
    C00095 -> R00866 -> C01094
    C00096 -> R00883 -> C00636
    C00097 -> R00895 -> C00957
    C00097 -> R00897 -> C00283
    C00097 -> R00897 -> C00979
    C00097 -> R00899 -> C01419
    C00097 -> R04859 -> C00979
    C00099 -> R00489 -> C00049
    C00099 -> R00907 -> C00222
    C00099 -> R02474 -> C00864
    C00100 -> R00919 -> C00894
    C00100 -> R00927 -> C03344
    C00100 -> R04432 -> C00894
    C00100 -> R10161 -> C00894
    C00101 -> R00937 -> C00504
    C00101 -> R00940 -> C00504
    C00101 -> R00943 -> C00234
    C00101 -> R01225 -> C00143
    C00101 -> R04326 -> C00445
    C00101 -> R06613 -> C00143
    C00103 -> R00287 -> C00029
    C00103 -> R00951 -> C00498
    C00103 -> R02111 -> C00369
    C00105 -> R00662 -> C00075
    C00105 -> R00964 -> C00299
    C00105 -> R00965 -> C01103
    C00105 -> R00967 -> C00299
    C00105 -> R00968 -> C00299
    C00105 -> R00970 -> C00299
    C00105 -> R01549 -> C00299
    C00105 -> R01880 -> C00299
    C00105 -> R02097 -> C00299
    C00105 -> R02327 -> C00299
    C00105 -> R02332 -> C00299
    C00106 -> R00974 -> C00380
    C00106 -> R00978 -> C00429
    C00108 -> R00985 -> C00251
    C00108 -> R00986 -> C00251
    C00109 -> R00996 -> C00188
    C00109 -> R00999 -> C01118
    C00111 -> R01011 -> C00184
    C00114 -> R01022 -> C00576
    C00114 -> R08557 -> C00576
    C00114 -> R08558 -> C00576
    C00116 -> R00841 -> C00093
    C00116 -> R01039 -> C00184
    C00118 -> R01063 -> C00236
    C00118 -> R01064 -> C01286
    C00120 -> R01078 -> C01909
    C00120 -> R10127 -> C20386
    C00122 -> R00490 -> C00049
    C00122 -> R01086 -> C03406
    C00123 -> R01088 -> C00233
    C00124 -> R01678 -> C00243
    C00129 -> R05884 -> C11811
    C00130 -> R01128 -> C00262
    C00131 -> R02014 -> C00002
    C00131 -> R11634 -> C00002
    C00134 -> R00670 -> C00077
    C00134 -> R01151 -> C00555
    C00134 -> R01156 -> C02714
    C00134 -> R01157 -> C00179
    C00134 -> R01920 -> C00170
    C00134 -> R08714 -> C00555
    C00135 -> R01163 -> C01929
    C00136 -> R01171 -> C00877
    C00136 -> R01176 -> C00246
    C00137 -> R01186 -> C03546
    C00137 -> R01187 -> C04006
    C00137 -> R07279 -> C04006
    C00140 -> R00022 -> C01674
    C00141 -> R01434 -> C00183
    C00141 -> R04441 -> C04272
    C00143 -> R01217 -> C00440
    C00143 -> R01218 -> C00445
    C00143 -> R01225 -> C00101
    C00143 -> R04125 -> C01242
    C00143 -> R06613 -> C00101
    C00144 -> R01230 -> C00655
    C00144 -> R01231 -> C00655
    C00147 -> R01245 -> C00212
    C00149 -> R00214 -> C00011
    C00149 -> R00216 -> C00011
    C00152 -> R00483 -> C00049
    C00154 -> R01278 -> C05272
    C00155 -> R01286 -> C02291
    C00155 -> R01287 -> C01077
    C00155 -> R01288 -> C01118
    C00155 -> R01290 -> C02291
    C00155 -> R01291 -> C03539
    C00155 -> R02026 -> C01077
    C00155 -> R02821 -> C00073
    C00155 -> R10305 -> C02291
    C00156 -> R05000 -> C05848
    C00157 -> R07064 -> C01595
    C00158 -> R00352 -> C00024
    C00158 -> R00352 -> C00036
    C00158 -> R01324 -> C00311
    C00159 -> R01326 -> C00275
    C00160 -> R00717 -> C00048
    C00160 -> R01334 -> C00988
    C00166 -> R00688 -> C00079
    C00166 -> R00689 -> C00079
    C00166 -> R00692 -> C00079
    C00169 -> R00149 -> C00014
    C00169 -> R00575 -> C00064
    C00170 -> R01920 -> C00134
    C00170 -> R01920 -> C01137
    C00173 -> R10121 -> C19845
    C00178 -> R01415 -> C21028
    C00178 -> R01570 -> C00214
    C00179 -> R00566 -> C00062
    C00181 -> R01433 -> C02352
    C00183 -> R01434 -> C00141
    C00184 -> R01039 -> C00116
    C00185 -> R02887 -> C01898
    C00185 -> R11308 -> C01898
    C00186 -> R00703 -> C00022
    C00186 -> R01447 -> C00022
    C00188 -> R01466 -> C01102
    C00191 -> R01478 -> C03033
    C00194 -> R01492 -> C00853
    C00194 -> R05223 -> C05775
    C00194 -> R05223 -> C06510
    C00198 -> R00305 -> C00031
    C00198 -> R06620 -> C00031
    C00198 -> R07147 -> C00031
    C00199 -> R01523 -> C01182
    C00199 -> R01528 -> C00345
    C00199 -> R10221 -> C00345
    C00204 -> R01540 -> C00817
    C00206 -> R02017 -> C00008
    C00208 -> R01555 -> C00031
    C00208 -> R11262 -> C01935
    C00212 -> R00183 -> C00020
    C00212 -> R00185 -> C00020
    C00214 -> R01569 -> C00364
    C00217 -> R01579 -> C00819
    C00222 -> R00907 -> C00099
    C00222 -> R09983 -> C20253
    C00224 -> R00508 -> C00053
    C00224 -> R00530 -> C00059
    C00224 -> R00531 -> C00059
    C00227 -> R00320 -> C00033
    C00232 -> R09281 -> C00989
    C00232 -> R10178 -> C00334
    C00233 -> R01088 -> C00123
    C00234 -> R00943 -> C00101
    C00235 -> R01122 -> C04432
    C00236 -> R01063 -> C00118
    C00239 -> R01667 -> C00705
    C00245 -> R01682 -> C00506
    C00246 -> R01176 -> C00136
    C00249 -> R01274 -> C00154
    C00250 -> R00173 -> C00018
    C00251 -> R01714 -> C01269
    C00253 -> R01268 -> C00153
    C00254 -> R01730 -> C01179
    C00254 -> R07276 -> C00826
    C00255 -> R00066 -> C04332
    C00255 -> R00548 -> C00061
    C00255 -> R00550 -> C00061
    C00255 -> R08574 -> C00061
    C00257 -> R01741 -> C06473
    C00262 -> R01128 -> C00130
    C00262 -> R01244 -> C00147
    C00262 -> R01769 -> C00385
    C00262 -> R02748 -> C05512
    C00267 -> R01678 -> C00243
    C00267 -> R01788 -> C00668
    C00267 -> R02189 -> C00668
    C00267 -> R09085 -> C00668
    C00269 -> R01799 -> C00416
    C00275 -> R01326 -> C00159
    C00275 -> R02630 -> C00159
    C00279 -> R01067 -> C00085
    C00279 -> R01829 -> C00447
    C00283 -> R00859 -> C00094
    C00283 -> R00861 -> C00094
    C00286 -> R01858 -> C00361
    C00286 -> R02020 -> C00044
    C00286 -> R11633 -> C00044
    C00294 -> R01126 -> C00130
    C00294 -> R01560 -> C00212
    C00295 -> R01867 -> C00337
    C00295 -> R01869 -> C00337
    C00299 -> R01878 -> C00475
    C00311 -> R00267 -> C00026
    C00311 -> R00709 -> C00026
    C00311 -> R01324 -> C00158
    C00314 -> R01911 -> C00627
    C00315 -> R01918 -> C05730
    C00315 -> R01920 -> C00134
    C00315 -> R01920 -> C01137
    C00327 -> R01398 -> C00077
    C00327 -> R01398 -> C00169
    C00327 -> R09107 -> C15532
    C00330 -> R01967 -> C00362
    C00334 -> R00261 -> C00025
    C00334 -> R01986 -> C00555
    C00334 -> R02549 -> C00555
    C00334 -> R07419 -> C15767
    C00334 -> R10178 -> C00232
    C00337 -> R01867 -> C00295
    C00337 -> R01869 -> C00295
    C00341 -> R01658 -> C00129
    C00341 -> R01658 -> C00235
    C00344 -> R02029 -> C03892
    C00345 -> R01737 -> C00257
    C00345 -> R02035 -> C01236
    C00350 -> R02055 -> C02737
    C00350 -> R07376 -> C02737
    C00352 -> R00768 -> C05345
    C00357 -> R01201 -> C00140
    C00360 -> R02089 -> C00559
    C00361 -> R02019 -> C00035
    C00362 -> R01967 -> C00330
    C00364 -> R01567 -> C00214
    C00364 -> R02101 -> C00365
    C00365 -> R02100 -> C00460
    C00366 -> R02103 -> C00385
    C00366 -> R02107 -> C00385
    C00369 -> R02110 -> C00718
    C00376 -> R02124 -> C00473
    C00376 -> R08379 -> C00473
    C00378 -> R02135 -> C01081
    C00380 -> R02137 -> C00475
    C00380 -> R02296 -> C00475
    C00385 -> R01676 -> C00242
    C00385 -> R01768 -> C00262
    C00385 -> R01769 -> C00262
    C00385 -> R02107 -> C00366
    C00387 -> R01227 -> C00144
    C00399 -> R08781 -> C17562
    C00407 -> R02197 -> C00671
    C00415 -> R02237 -> C00921
    C00416 -> R02239 -> C00641
    C00416 -> R02240 -> C00641
    C00416 -> R02241 -> C00681
    C00416 -> R09944 -> C00641
    C00429 -> R00978 -> C00106
    C00430 -> R02272 -> C03741
    C00437 -> R02282 -> C00077
    C00437 -> R02283 -> C01250
    C00438 -> R01397 -> C00169
    C00440 -> R01217 -> C00143
    C00440 -> R01224 -> C00143
    C00440 -> R07168 -> C00143
    C00441 -> R02291 -> C03082
    C00445 -> R01218 -> C00143
    C00445 -> R02301 -> C03479
    C00446 -> R00502 -> C00052
    C00446 -> R01092 -> C00984
    C00447 -> R01829 -> C00279
    C00448 -> R02003 -> C00129
    C00448 -> R02003 -> C00341
    C00455 -> R00103 -> C00003
    C00455 -> R02324 -> C03150
    C00458 -> R02022 -> C00063
    C00458 -> R11636 -> C00063
    C00460 -> R02325 -> C00458
    C00469 -> R00746 -> C00084
    C00469 -> R05198 -> C00084
    C00469 -> R09127 -> C00084
    C00470 -> R02362 -> C00714
    C00473 -> R02124 -> C00376
    C00473 -> R08379 -> C00376
    C00475 -> R02296 -> C00380
    C00493 -> R02415 -> C02637
    C00498 -> R00951 -> C00103
    C00519 -> R02466 -> C00606
    C00522 -> R02474 -> C00864
    C00534 -> R02494 -> C00647
    C00541 -> R00097 -> C00992
    C00541 -> R00107 -> C00992
    C00555 -> R01151 -> C00134
    C00555 -> R01155 -> C00134
    C00555 -> R01986 -> C00334
    C00555 -> R08714 -> C00134
    C00558 -> R01983 -> C00333
    C00559 -> R02089 -> C00360
    C00576 -> R01022 -> C00114
    C00576 -> R08211 -> C00719
    C00576 -> R08557 -> C00114
    C00576 -> R08558 -> C00114
    C00590 -> R02596 -> C15805
    C00601 -> R02537 -> C07086
    C00624 -> R00259 -> C00025
    C00627 -> R01909 -> C00314
    C00627 -> R01911 -> C00314
    C00627 -> R05838 -> C11638
    C00631 -> R08572 -> C00258
    C00636 -> R00883 -> C00096
    C00641 -> R02239 -> C00416
    C00641 -> R09944 -> C00416
    C00645 -> R00414 -> C00043
    C00647 -> R02493 -> C00534
    C00647 -> R02494 -> C00534
    C00655 -> R01130 -> C00130
    C00666 -> R02734 -> C04421
    C00668 -> R01788 -> C00267
    C00668 -> R02189 -> C00267
    C00668 -> R09085 -> C00267
    C00669 -> R00894 -> C00025
    C00669 -> R00894 -> C00097
    C00671 -> R02197 -> C00407
    C00671 -> R05070 -> C06007
    C00681 -> R00851 -> C00093
    C00689 -> R02737 -> C00029
    C00692 -> R02783 -> C01212
    C00705 -> R01667 -> C00239
    C00705 -> R02024 -> C00112
    C00718 -> R02421 -> C00498
    C00719 -> R02565 -> C00576
    C00719 -> R02566 -> C00576
    C00719 -> R08211 -> C00576
    C00721 -> R02112 -> C00369
    C00748 -> R02864 -> C05778
    C00760 -> R02889 -> C00029
    C00798 -> R01908 -> C00313
    C00817 -> R02555 -> C00558
    C00826 -> R01731 -> C00254
    C00826 -> R07276 -> C00254
    C00828 -> R04993 -> C05818
    C00857 -> R00257 -> C00003
    C00860 -> R03013 -> C01100
    C00864 -> R02473 -> C00099
    C00864 -> R02473 -> C00522
    C00864 -> R02474 -> C00099
    C00864 -> R02474 -> C00522
    C00877 -> R01171 -> C00136
    C00877 -> R01175 -> C00136
    C00894 -> R04432 -> C00100
    C00894 -> R10161 -> C00100
    C00921 -> R03066 -> C00568
    C00921 -> R03066 -> C01300
    C00921 -> R03067 -> C00568
    C00921 -> R03067 -> C04807
    C00931 -> R00036 -> C00430
    C00944 -> R03083 -> C04691
    C00957 -> R00895 -> C00097
    C00966 -> R01226 -> C00141
    C00979 -> R04859 -> C00097
    C00984 -> R10619 -> C00124
    C00989 -> R09281 -> C00232
    C00992 -> R00107 -> C00541
    C00993 -> R01150 -> C00133
    C01007 -> R05707 -> C00255
    C01007 -> R09750 -> C00255
    C01013 -> R09289 -> C00222
    C01024 -> R00084 -> C00931
    C01037 -> R03231 -> C01092
    C01037 -> R10699 -> C01092
    C01050 -> R03191 -> C04631
    C01050 -> R03192 -> C04631
    C01051 -> R03165 -> C01024
    C01077 -> R01287 -> C00155
    C01079 -> R03220 -> C03263
    C01079 -> R03222 -> C02191
    C01079 -> R06895 -> C03263
    C01081 -> R00615 -> C00068
    C01081 -> R02134 -> C00378
    C01081 -> R03223 -> C04327
    C01081 -> R03223 -> C04752
    C01081 -> R10712 -> C04752
    C01081 -> R10712 -> C20247
    C01083 -> R02727 -> C00031
    C01083 -> R02778 -> C00689
    C01092 -> R03210 -> C01063
    C01092 -> R10124 -> C19845
    C01092 -> R10699 -> C01037
    C01094 -> R00866 -> C00095
    C01094 -> R03232 -> C00095
    C01097 -> R05571 -> C06311
    C01102 -> R01771 -> C00263
    C01134 -> R03269 -> C04352
    C01137 -> R00178 -> C00019
    C01146 -> R00013 -> C00048
    C01146 -> R01394 -> C00168
    C01157 -> R03291 -> C04281
    C01157 -> R03293 -> C04281
    C01157 -> R03295 -> C04281
    C01165 -> R03313 -> C03287
    C01170 -> R00420 -> C00043
    C01179 -> R00729 -> C00082
    C01179 -> R01730 -> C00254
    C01179 -> R09254 -> C00082
    C01180 -> R00648 -> C00073
    C01182 -> R01523 -> C00199
    C01185 -> R03004 -> C00857
    C01185 -> R03348 -> C03722
    C01187 -> R03350 -> C04478
    C01190 -> R03354 -> C01290
    C01190 -> R03355 -> C01290
    C01212 -> R03193 -> C01050
    C01216 -> R03033 -> C00880
    C01222 -> R00888 -> C00096
    C01236 -> R02736 -> C01172
    C01236 -> R10907 -> C01172
    C01242 -> R03425 -> C00037
    C01250 -> R03443 -> C04133
    C01267 -> R03457 -> C04666
    C01268 -> R03459 -> C01304
    C01286 -> R03387 -> C01216
    C01290 -> R03354 -> C01190
    C01300 -> R03504 -> C04874
    C01304 -> R00425 -> C00044
    C01419 -> R00494 -> C00051
    C01528 -> R03599 -> C05688
    C01595 -> R07064 -> C00157
    C01672 -> R00462 -> C00047
    C01672 -> R06740 -> C12455
    C01674 -> R02334 -> C00461
    C01762 -> R02719 -> C00655
    C01832 -> R03856 -> C03221
    C01847 -> R05705 -> C00061
    C01847 -> R05706 -> C00061
    C01898 -> R11308 -> C00185
    C01909 -> R03182 -> C01037
    C01929 -> R03012 -> C00860
    C01944 -> R03776 -> C05276
    C02059 -> R06859 -> C13309
    C02191 -> R03222 -> C01079
    C02191 -> R09489 -> C01079
    C02222 -> R08121 -> C16476
    C02232 -> R06941 -> C14145
    C02282 -> R03652 -> C00064
    C02291 -> R01290 -> C00155
    C02291 -> R03217 -> C01077
    C02291 -> R03260 -> C01118
    C02291 -> R10305 -> C00155
    C02325 -> R03919 -> C15806
    C02463 -> R03194 -> C01051
    C02504 -> R01213 -> C00141
    C02565 -> R02922 -> C00791
    C02593 -> R03989 -> C05273
    C02637 -> R02415 -> C00493
    C02646 -> R04007 -> C15804
    C02714 -> R01154 -> C00134
    C02714 -> R01156 -> C00134
    C02723 -> R04027 -> C06178
    C02730 -> R04031 -> C05817
    C02737 -> R01800 -> C00269
    C02737 -> R07376 -> C00350
    C02741 -> R04035 -> C02739
    C02987 -> R05578 -> C00025
    C03028 -> R00616 -> C00068
    C03028 -> R00618 -> C00068
    C03028 -> R11319 -> C00068
    C03082 -> R00480 -> C00049
    C03089 -> R01401 -> C00170
    C03090 -> R01072 -> C00119
    C03160 -> R04030 -> C02730
    C03175 -> R02412 -> C00493
    C03221 -> R03856 -> C01832
    C03221 -> R03857 -> C01832
    C03263 -> R03197 -> C01051
    C03287 -> R00239 -> C00025
    C03291 -> R07125 -> C14899
    C03296 -> R00832 -> C00062
    C03373 -> R04208 -> C04640
    C03406 -> R01954 -> C00049
    C03406 -> R01954 -> C00327
    C03415 -> R04189 -> C03296
    C03453 -> R03966 -> C02501
    C03492 -> R03018 -> C00864
    C03492 -> R04230 -> C04352
    C03539 -> R00194 -> C00021
    C03589 -> R02601 -> C00596
    C03657 -> R07262 -> C15547
    C03684 -> R04286 -> C04895
    C03722 -> R04292 -> C05840
    C03741 -> R04109 -> C02987
    C03794 -> R01135 -> C00130
    C03838 -> R04144 -> C03090
    C03892 -> R01801 -> C00269
    C03912 -> R01253 -> C00148
    C03972 -> R04198 -> C20258
    C03972 -> R04199 -> C20258
    C04006 -> R07279 -> C00137
    C04043 -> R04300 -> C03758
    C04121 -> R03351 -> C01187
    C04133 -> R02649 -> C00624
    C04257 -> R02705 -> C00645
    C04281 -> R03295 -> C01157
    C04317 -> R03438 -> C05212
    C04317 -> R07387 -> C05212
    C04327 -> R04448 -> C04294
    C04332 -> R04457 -> C04732
    C04332 -> R04457 -> C15556
    C04352 -> R04230 -> C03492
    C04352 -> R04231 -> C03492
    C04376 -> R04325 -> C03838
    C04405 -> R04204 -> C03345
    C04432 -> R01122 -> C00235
    C04442 -> R01541 -> C00204
    C04442 -> R02036 -> C00345
    C04454 -> R03458 -> C01268
    C04462 -> R04365 -> C03972
    C04478 -> R03254 -> C01112
    C04501 -> R05332 -> C06156
    C04556 -> R03471 -> C01279
    C04556 -> R03472 -> C03373
    C04631 -> R00660 -> C00043
    C04635 -> R04413 -> C04756
    C04635 -> R07379 -> C04756
    C04640 -> R04463 -> C04376
    C04652 -> R04550 -> C06022
    C04666 -> R04558 -> C04916
    C04677 -> R06975 -> C04734
    C04691 -> R01826 -> C00074
    C04691 -> R01826 -> C00279
    C04702 -> R04573 -> C05892
    C04732 -> R07280 -> C04454
    C04734 -> R04560 -> C04677
    C04734 -> R06975 -> C04677
    C04738 -> R04567 -> C00043
    C04752 -> R04509 -> C04556
    C04756 -> R04413 -> C04635
    C04778 -> R04148 -> C03114
    C04807 -> R03503 -> C01300
    C04824 -> R04549 -> C04652
    C04851 -> R05629 -> C04702
    C04874 -> R04620 -> C04895
    C04895 -> R04639 -> C06148
    C04896 -> R04037 -> C02741
    C04916 -> R04640 -> C04896
    C04919 -> R04657 -> C04932
    C04932 -> R04606 -> C04652
    C04932 -> R04606 -> C04824
    C05172 -> R03595 -> C01528
    C05212 -> R03438 -> C04317
    C05270 -> R06985 -> C05271
    C05271 -> R04751 -> C05270
    C05271 -> R06985 -> C05270
    C05272 -> R01278 -> C00154
    C05272 -> R01279 -> C00154
    C05273 -> R03989 -> C02593
    C05273 -> R03990 -> C02593
    C05274 -> R04753 -> C05275
    C05275 -> R04753 -> C05274
    C05275 -> R04754 -> C05274
    C05276 -> R03776 -> C01944
    C05276 -> R03777 -> C01944
    C05345 -> R00867 -> C00095
    C05345 -> R02073 -> C05378
    C05345 -> R04780 -> C05378
    C05345 -> R09084 -> C05378
    C05378 -> R02073 -> C05345
    C05378 -> R04779 -> C05345
    C05378 -> R09084 -> C05345
    C05527 -> R02619 -> C00606
    C05730 -> R01917 -> C00051
    C05730 -> R01917 -> C00315
    C05744 -> R04355 -> C01209
    C05744 -> R04355 -> C03939
    C05744 -> R10707 -> C01209
    C05746 -> R04952 -> C01209
    C05746 -> R04952 -> C05745
    C05750 -> R04957 -> C01209
    C05750 -> R04957 -> C05749
    C05753 -> R04960 -> C01209
    C05753 -> R04960 -> C05752
    C05756 -> R04963 -> C01209
    C05756 -> R04963 -> C05755
    C05759 -> R04726 -> C01209
    C05759 -> R04726 -> C05223
    C05762 -> R04968 -> C01209
    C05762 -> R04968 -> C05761
    C05768 -> R04972 -> C05766
    C05775 -> R04594 -> C04778
    C05778 -> R03947 -> C02463
    C05791 -> R04979 -> C05787
    C05807 -> R04985 -> C05848
    C05807 -> R08768 -> C17551
    C05817 -> R08166 -> C16519
    C05818 -> R05617 -> C03657
    C05840 -> R00481 -> C00049
    C05840 -> R07407 -> C00049
    C05840 -> R07410 -> C00049
    C05848 -> R04985 -> C05807
    C05848 -> R05000 -> C00156
    C05893 -> R05662 -> C04851
    C05921 -> R01074 -> C00120
    C05922 -> R00428 -> C00044
    C05923 -> R05046 -> C05922
    C05924 -> R09395 -> C18239
    C05931 -> R05049 -> C05932
    C05932 -> R04217 -> C03415
    C05946 -> R05052 -> C05947
    C05946 -> R05053 -> C05947
    C05947 -> R04444 -> C04281
    C05947 -> R04445 -> C04281
    C05947 -> R05053 -> C05946
    C05980 -> R07390 -> C00344
    C05980 -> R11062 -> C00344
    C06000 -> R04224 -> C03460
    C06001 -> R05066 -> C06002
    C06002 -> R05066 -> C06001
    C06006 -> R08648 -> C00022
    C06006 -> R08648 -> C00109
    C06010 -> R00226 -> C00022
    C06022 -> R04587 -> C04738
    C06024 -> R04658 -> C04121
    C06024 -> R04658 -> C04919
    C06025 -> R05074 -> C04121
    C06025 -> R05074 -> C06024
    C06026 -> R05075 -> C06251
    C06040 -> R05081 -> C06041
    C06041 -> R05081 -> C06040
    C06148 -> R05048 -> C05923
    C06157 -> R01940 -> C00322
    C06178 -> R04027 -> C02723
    C06250 -> R05145 -> C05921
    C06251 -> R05146 -> C06025
    C06311 -> R05570 -> C01697
    C06329 -> R05511 -> C04706
    C06397 -> R05644 -> C07838
    C06398 -> R05176 -> C06397
    C06427 -> R07859 -> C00157
    C06427 -> R07860 -> C00157
    C06473 -> R01741 -> C00257
    C06506 -> R05220 -> C06505
    C06509 -> R05221 -> C06508
    C06509 -> R06558 -> C06508
    C06510 -> R05222 -> C06509
    C06613 -> R05233 -> C06611
    C07086 -> R02536 -> C00601
    C07086 -> R02537 -> C00601
    C07335 -> R05681 -> C06055
    C07479 -> R05389 -> C07478
    C07836 -> R05645 -> C05382
    C07838 -> R05647 -> C11472
    C11434 -> R05688 -> C11437
    C11435 -> R05633 -> C11434
    C11436 -> R05634 -> C11435
    C11453 -> R05637 -> C11436
    C11472 -> R05646 -> C07836
    C11811 -> R08689 -> C11453
    C11811 -> R10859 -> C11453
    C12248 -> R06601 -> C11821
    C12455 -> R06740 -> C01672
    C12835 -> R06835 -> C12834
    C13309 -> R06858 -> C03657
    C14145 -> R06942 -> C14144
    C14899 -> R07677 -> C16186
    C15556 -> R07281 -> C00199
    C15667 -> R07404 -> C03373
    C15672 -> R07411 -> C00032
    C15699 -> R07414 -> C00134
    C15700 -> R07415 -> C15699
    C15767 -> R07417 -> C15700
    C15767 -> R07418 -> C15700
    C15804 -> R04007 -> C02646
    C15805 -> R02596 -> C00590
    C15806 -> R03919 -> C02325
    C15809 -> R10246 -> C00082
    C15812 -> R07460 -> C00097
    C15812 -> R07460 -> C15811
    C15813 -> R07459 -> C15810
    C15814 -> R07461 -> C15812
    C15814 -> R07461 -> C15813
    C15996 -> R09978 -> C20248
    C16186 -> R07671 -> C00072
    C16236 -> R07766 -> C05752
    C16237 -> R07767 -> C16236
    C16237 -> R07769 -> C16239
    C16237 -> R11143 -> C16241
    C16239 -> R07768 -> C05752
    C16331 -> R07891 -> C16330
    C16335 -> R07895 -> C16334
    C16339 -> R07899 -> C16338
    C16348 -> R05234 -> C06612
    C16519 -> R08165 -> C00885
    C16675 -> R07605 -> C15996
    C17551 -> R08768 -> C05807
    C17551 -> R08769 -> C17552
    C17552 -> R08769 -> C17551
    C17552 -> R08773 -> C17560
    C17560 -> R08773 -> C17552
    C17560 -> R08774 -> C17561
    C17561 -> R08774 -> C17560
    C17561 -> R08775 -> C17562
    C17562 -> R08775 -> C17561
    C17562 -> R08781 -> C00399
    C18237 -> R09735 -> C19848
    C18239 -> R11372 -> C21310
    C19673 -> R09543 -> C01209
    C19845 -> R09725 -> C19846
    C19845 -> R10121 -> C00173
    C19846 -> R10122 -> C20378
    C19848 -> R09726 -> C05924
    C19871 -> R11581 -> C18237
    C20231 -> R09936 -> C00106
    C20239 -> R09959 -> C04895
    C20246 -> R10247 -> C11437
    C20246 -> R10247 -> C15809
    C20246 -> R10247 -> C15814
    C20248 -> R10002 -> C20239
    C20249 -> R09947 -> C20231
    C20253 -> R09982 -> C20249
    C20258 -> R10147 -> C00441
    C20372 -> R10115 -> C01209
    C20372 -> R10115 -> C19673
    C20373 -> R10116 -> C20372
    C20374 -> R10117 -> C20373
    C20375 -> R10118 -> C20374
    C20376 -> R10119 -> C01209
    C20376 -> R10119 -> C20375
    C20377 -> R10120 -> C20376
    C21015 -> R10993 -> C02356
    C21016 -> R10994 -> C21015
    C21028 -> R01415 -> C00178
    C21284 -> R11329 -> C05770
    C21310 -> R09394 -> C00044
    G09660 -> R07818 -> G13040
    G13040 -> R07818 -> G09660

Conclusion
----------
Many reactions occur in the overview map eco01100, although they exist no where in the union of all known metabolism pathways.

Further research shows the reasons for such strange results:

- If a gene of eco produces an enzyme that can be used in multiple pathways, 
  eco01100 automatically includes the reactions performed by that enzyme into all possible positions of the global map, 
  even those positions inside a pathway which eco supposedly does not have. 
  Because the reaction acts on a different substrate in this "non-existing" pathway, it stands out in this experiment.
  This is most likely not an error in KEGG's data, but a corner-case where our methodology fails to handle the data.
  Example: G09660 -> R07818 -> G13040 (gene eco:b1617)
  
- The same reaction (same substrate and product) is marked as irreversible in a single pathway,
  in eco01100, however, it is marked a reversible. 01100 generally seems to not support irreversible edges.
  Because this experiment uses a directed graph, the reverse direction stands out.
  Example: C21310 -> R09394 -> C00044
  
- Maybe more.

The overview map 01100 is, therefore, not suitable for exact research of reactions.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco01100 = eco.getPathway('01100')
    allNonOverviewPathways = eco.getMetabolicPathways(includeOverviewMaps = False)
    
    #- Convert to substance-reaction graph.
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    allNonOverviewPathways_reactionGraph = SubstanceReactionGraph.fromPathway(allNonOverviewPathways)
    
    #- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
    difference_reactionGraph = eco01100_reactionGraph.difference(allNonOverviewPathways_reactionGraph, subtractNodes = False)
    
    #- Print reactions (substrate -> reaction ID -> product).
    output = []
    edges = difference_reactionGraph.getEdges()
    for edge in edges:
        substrate, product, reactionID = edge
        output.append(substrate.__str__() + ' -> ' + reactionID.__str__() + ' -> ' + product.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
