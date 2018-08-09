"""
Question
--------
Assuming eco and sso had a common ancestor, and further assuming eco has evolved into a more versatile metabolism.
Which functions - defined by a new EC number - occuring in eco, but not in sso, may have originated from a neofunctionalisation? 

Method
------
Similar to :mod:`17`:
- Get all known pathways of eco and sso from KEGG.
- Convert each pathway into a substance-ecNumber graph. Using a substance-enzyme graph as intermediate. This is the (a little incomplete) metabolic network.
- Remove multifunctional enzymes, meaning enzymes associated with more than one EC number. Helps to reduce false gene duplications.
- Combine both species' networks to a consensus network, by INTERSECTION operation.
- Create the networks of EC numbers that only occur in eco.
Similar to :mod:`17`, without triangle condition:
- For each eco-only EC number, gather any eco gene which encodes an enzyme associated with this EC number.
- For each such gene, find homologs within eco (paralogs, threshold=200).

- If there exists at least one paralog with a different function, report a neofunctionalisation, in the following format:
  paralog with different EC number (EC number) + neofunctionalised gene (new function) 
  This may report different paralog groups per neofunctionalised gene, when the group of paralogs contains multiple functions. May also report other neofunctionalised genes as paralogs, because a single group of paralogs can have multiple neofunctionalisations.

- Report percentage neofunctionalisations of all new functions.

Result
------

::

    323 results
    eco:b0063 (2.7.1.16) + eco:b3564 (2.7.1.17)
    eco:b0077 (2.2.1.6) + eco:b0507 (4.1.1.47)
    eco:b0077 (2.2.1.6) + eco:b2373 (4.1.1.8)
    eco:b0085 (6.3.2.13) + eco:b0086 (6.3.2.10)
    eco:b0086 (6.3.2.10) + eco:b0085 (6.3.2.13)
    eco:b0086 (6.3.2.10) + eco:b0088 (6.3.2.9)
    eco:b0086 (6.3.2.10) + eco:b0091 (6.3.2.8)
    eco:b0088 (6.3.2.9) + eco:b0086 (6.3.2.10)
    eco:b0088 (6.3.2.9) + eco:b0091 (6.3.2.8)
    eco:b0091 (6.3.2.8) + eco:b0086 (6.3.2.10)
    eco:b0091 (6.3.2.8) + eco:b0088 (6.3.2.9)
    eco:b0104 (1.7.1.7) + eco:b2508 (1.1.1.205)
    eco:b0115 (2.3.1.12) + eco:b0727 (2.3.1.61)
    eco:b0116 (1.8.1.4) + eco:b3365 (1.7.1.15)
    eco:b0116 (1.8.1.4) + eco:b3500 (1.8.1.7)
    eco:b0116 (1.8.1.4) + eco:b3962 (1.6.1.1)
    eco:b0154 (5.4.3.8) + eco:b0774 (2.6.1.62)
    eco:b0154 (5.4.3.8) + eco:b1748 (2.6.1.81)
    eco:b0154 (5.4.3.8) + eco:b3073 (2.6.1.82)
    eco:b0179 (2.3.1.191) + eco:b0181 (2.3.1.129)
    eco:b0181 (2.3.1.129) + eco:b0179 (2.3.1.191)
    eco:b0186 (4.1.1.18) + eco:b0693 (4.1.1.17)
    eco:b0186 (4.1.1.18) + eco:b2965 (4.1.1.17)
    eco:b0186 (4.1.1.18) + eco:b4117 (4.1.1.19)
    eco:b0261 (2.1.1.10) + eco:b4019 (2.1.1.13)
    eco:b0268 (4.3.3.7) + eco:b3225 (4.1.3.3)
    eco:b0312 (1.2.1.8) + eco:b1300 (1.2.1.99)
    eco:b0312 (1.2.1.8) + eco:b1385 (1.2.1.39)
    eco:b0312 (1.2.1.8) + eco:b1444 (1.2.1.19)
    eco:b0312 (1.2.1.8) + eco:b1746 (1.2.1.71)
    eco:b0312 (1.2.1.8) + eco:b3588 (1.2.1.-)
    eco:b0335 (6.2.1.17) + eco:b0586 (6.3.2.14)
    eco:b0394 (2.7.1.4) + eco:b1119 (2.7.1.59)
    eco:b0394 (2.7.1.4) + eco:b3222 (2.7.1.60)
    eco:b0394 (2.7.1.4) + eco:b4084 (2.7.1.55)
    eco:b0403 (3.2.1.20) + eco:b4239 (3.2.1.93)
    eco:b0507 (4.1.1.47) + eco:b2373 (4.1.1.8)
    eco:b0512 (3.5.2.5) + eco:b2873 (3.5.2.2)
    eco:b0517 (1.1.1.350) + eco:b3575 (1.1.1.130)
    eco:b0586 (6.3.2.14) + eco:b0335 (6.2.1.17)
    eco:b0586 (6.3.2.14) + eco:b2260 (6.2.1.26)
    eco:b0593 (5.4.4.2) + eco:b1812 (2.6.1.85)
    eco:b0596 (1.3.1.28) + eco:b2842 (1.1.1.127)
    eco:b0678 (3.5.99.6) + eco:b3141 (5.3.1.-)
    eco:b0679 (2.7.1.193) + eco:b1101 (2.7.1.199)
    eco:b0679 (2.7.1.193) + eco:b2417 (2.7.1.-)
    eco:b0688 (5.4.2.2) + eco:b3176 (5.4.2.10)
    eco:b0693 (4.1.1.17) + eco:b0186 (4.1.1.18)
    eco:b0693 (4.1.1.17) + eco:b4117 (4.1.1.19)
    eco:b0693 (4.1.1.17) + eco:b4131 (4.1.1.18)
    eco:b0720 (2.3.3.1) + eco:b0333 (2.3.3.5)
    eco:b0759 (5.1.3.2) + eco:b2053 (4.2.1.47)
    eco:b0759 (5.1.3.2) + eco:b3619 (5.1.3.20)
    eco:b0774 (2.6.1.62) + eco:b1748 (2.6.1.81)
    eco:b0774 (2.6.1.62) + eco:b3073 (2.6.1.82)
    eco:b0776 (2.3.1.47) + eco:b3617 (2.3.1.29)
    eco:b0871 (1.2.5.1) + eco:b0507 (4.1.1.47)
    eco:b0871 (1.2.5.1) + eco:b2373 (4.1.1.8)
    eco:b0894 (1.8.5.3) + eco:b1587 (1.97.1.9)
    eco:b0894 (1.8.5.3) + eco:b1588 (1.97.1.9)
    eco:b0894 (1.8.5.3) + eco:b2206 (1.7.99.-)
    eco:b0894 (1.8.5.3) + eco:b3551 (1.-.-.-)
    eco:b0894 (1.8.5.3) + eco:b3894 (1.2.1.2)
    eco:b0908 (2.5.1.19) + eco:b3189 (2.5.1.7)
    eco:b0928 (2.6.1.1) + eco:b4054 (2.6.1.57)
    eco:b0935 (1.14.14.5) + eco:b1012 (1.14.99.46)
    eco:b1012 (1.14.99.46) + eco:b0935 (1.14.14.5)
    eco:b1054 (2.3.1.241) + eco:b1855 (2.3.1.243)
    eco:b1093 (1.1.1.100) + eco:b0596 (1.3.1.28)
    eco:b1093 (1.1.1.100) + eco:b2705 (1.1.1.140)
    eco:b1093 (1.1.1.100) + eco:b2842 (1.1.1.127)
    eco:b1095 (2.3.1.179) + eco:b2323 (2.3.1.41)
    eco:b1101 (2.7.1.199) + eco:b0679 (2.7.1.193)
    eco:b1119 (2.7.1.59) + eco:b0394 (2.7.1.4)
    eco:b1119 (2.7.1.59) + eco:b3222 (2.7.1.60)
    eco:b1119 (2.7.1.59) + eco:b4084 (2.7.1.55)
    eco:b1264 (4.1.3.27) + eco:b0593 (5.4.4.2)
    eco:b1264 (4.1.3.27) + eco:b1812 (2.6.1.85)
    eco:b1264 (4.1.3.27) + eco:b2265 (5.4.4.2)
    eco:b1300 (1.2.1.99) + eco:b0312 (1.2.1.8)
    eco:b1300 (1.2.1.99) + eco:b1385 (1.2.1.39)
    eco:b1300 (1.2.1.99) + eco:b1444 (1.2.1.19)
    eco:b1300 (1.2.1.99) + eco:b1746 (1.2.1.71)
    eco:b1300 (1.2.1.99) + eco:b3588 (1.2.1.-)
    eco:b1302 (2.6.1.19) + eco:b0774 (2.6.1.62)
    eco:b1302 (2.6.1.19) + eco:b1748 (2.6.1.81)
    eco:b1302 (2.6.1.19) + eco:b3073 (2.6.1.82)
    eco:b1309 (2.4.1.7) + eco:b4239 (3.2.1.93)
    eco:b1380 (1.1.1.28) + eco:b2320 (1.1.1.290)
    eco:b1385 (1.2.1.39) + eco:b0312 (1.2.1.8)
    eco:b1385 (1.2.1.39) + eco:b1300 (1.2.1.99)
    eco:b1385 (1.2.1.39) + eco:b1444 (1.2.1.19)
    eco:b1385 (1.2.1.39) + eco:b1746 (1.2.1.71)
    eco:b1385 (1.2.1.39) + eco:b3588 (1.2.1.-)
    eco:b1393 (4.2.1.17) + eco:b1394 (5.3.3.18)
    eco:b1393 (4.2.1.17) + eco:b2262 (4.1.3.36)
    eco:b1393 (4.2.1.17) + eco:b2919 (4.1.1.41)
    eco:b1394 (5.3.3.18) + eco:b1393 (4.2.1.17)
    eco:b1394 (5.3.3.18) + eco:b2262 (4.1.3.36)
    eco:b1394 (5.3.3.18) + eco:b2919 (4.1.1.41)
    eco:b1406 (1.1.1.65) + eco:b1771 (1.1.1.-)
    eco:b1444 (1.2.1.19) + eco:b0312 (1.2.1.8)
    eco:b1444 (1.2.1.19) + eco:b1300 (1.2.1.99)
    eco:b1444 (1.2.1.19) + eco:b1385 (1.2.1.39)
    eco:b1444 (1.2.1.19) + eco:b1746 (1.2.1.71)
    eco:b1444 (1.2.1.19) + eco:b3588 (1.2.1.-)
    eco:b1478 (1.1.1.1) + eco:b1580 (1.1.1.380)
    eco:b1478 (1.1.1.1) + eco:b2091 (1.1.1.251)
    eco:b1478 (1.1.1.1) + eco:b3616 (1.1.1.103)
    eco:b1479 (1.1.1.38) + eco:b2463 (1.1.1.40)
    eco:b1521 (1.1.1.58) + eco:b4323 (1.1.1.57)
    eco:b1580 (1.1.1.380) + eco:b2091 (1.1.1.251)
    eco:b1580 (1.1.1.380) + eco:b3616 (1.1.1.103)
    eco:b1581 (4.2.1.8) + eco:b2247 (4.2.1.90)
    eco:b1587 (1.97.1.9) + eco:b0894 (1.8.5.3)
    eco:b1587 (1.97.1.9) + eco:b2206 (1.7.99.-)
    eco:b1587 (1.97.1.9) + eco:b3551 (1.-.-.-)
    eco:b1587 (1.97.1.9) + eco:b3894 (1.2.1.2)
    eco:b1588 (1.97.1.9) + eco:b0894 (1.8.5.3)
    eco:b1588 (1.97.1.9) + eco:b2206 (1.7.99.-)
    eco:b1588 (1.97.1.9) + eco:b3551 (1.-.-.-)
    eco:b1588 (1.97.1.9) + eco:b3894 (1.2.1.2)
    eco:b1611 (4.2.1.2) + eco:b4139 (4.3.1.1)
    eco:b1612 (4.2.1.2) + eco:b3061 (4.2.1.32)
    eco:b1612 (4.2.1.2) + eco:b3062 (4.2.1.32)
    eco:b1617 (3.2.1.31) + eco:b0344 (3.2.1.23)
    eco:b1617 (3.2.1.31) + eco:b3076 (3.2.1.23)
    eco:b1723 (2.7.1.11) + eco:b2168 (2.7.1.56)
    eco:b1734 (3.2.1.86) + eco:b4119 (3.2.1.22)
    eco:b1746 (1.2.1.71) + eco:b0312 (1.2.1.8)
    eco:b1746 (1.2.1.71) + eco:b1300 (1.2.1.99)
    eco:b1746 (1.2.1.71) + eco:b1385 (1.2.1.39)
    eco:b1746 (1.2.1.71) + eco:b1444 (1.2.1.19)
    eco:b1746 (1.2.1.71) + eco:b3588 (1.2.1.-)
    eco:b1748 (2.6.1.81) + eco:b0774 (2.6.1.62)
    eco:b1748 (2.6.1.81) + eco:b3073 (2.6.1.82)
    eco:b1771 (1.1.1.-) + eco:b1406 (1.1.1.65)
    eco:b1779 (1.2.1.12) + eco:b2927 (1.2.1.72)
    eco:b1805 (6.2.1.3) + eco:b0335 (6.2.1.17)
    eco:b1805 (6.2.1.3) + eco:b0586 (6.3.2.14)
    eco:b1805 (6.2.1.3) + eco:b2260 (6.2.1.26)
    eco:b1812 (2.6.1.85) + eco:b0593 (5.4.4.2)
    eco:b1812 (2.6.1.85) + eco:b2265 (5.4.4.2)
    eco:b1817 (2.7.1.191) + eco:b3133 (2.7.1.-)
    eco:b1817 (2.7.1.191) + eco:b3138 (2.7.1.-)
    eco:b1855 (2.3.1.243) + eco:b1054 (2.3.1.241)
    eco:b2028 (1.1.1.22) + eco:b3787 (1.1.1.336)
    eco:b2041 (4.2.1.46) + eco:b2052 (1.1.1.271)
    eco:b2041 (4.2.1.46) + eco:b2053 (4.2.1.47)
    eco:b2048 (5.4.2.8) + eco:b3176 (5.4.2.10)
    eco:b2091 (1.1.1.251) + eco:b1580 (1.1.1.380)
    eco:b2091 (1.1.1.251) + eco:b3616 (1.1.1.103)
    eco:b2146 (1.3.1.1) + eco:b2878 (1.97.1.9)
    eco:b2166 (2.7.1.83) + eco:b2168 (2.7.1.56)
    eco:b2168 (2.7.1.56) + eco:b1723 (2.7.1.11)
    eco:b2168 (2.7.1.56) + eco:b2166 (2.7.1.83)
    eco:b2169 (2.7.1.202) + eco:b3599 (2.7.1.197)
    eco:b2206 (1.7.99.-) + eco:b0894 (1.8.5.3)
    eco:b2206 (1.7.99.-) + eco:b1587 (1.97.1.9)
    eco:b2206 (1.7.99.-) + eco:b1588 (1.97.1.9)
    eco:b2206 (1.7.99.-) + eco:b3551 (1.-.-.-)
    eco:b2206 (1.7.99.-) + eco:b3894 (1.2.1.2)
    eco:b2245 (4.1.2.53) + eco:b3126 (4.1.2.20)
    eco:b2247 (4.2.1.90) + eco:b1581 (4.2.1.8)
    eco:b2260 (6.2.1.26) + eco:b0586 (6.3.2.14)
    eco:b2262 (4.1.3.36) + eco:b1393 (4.2.1.17)
    eco:b2262 (4.1.3.36) + eco:b1394 (5.3.3.18)
    eco:b2262 (4.1.3.36) + eco:b2919 (4.1.1.41)
    eco:b2265 (5.4.4.2) + eco:b1812 (2.6.1.85)
    eco:b2296 (2.7.2.1) + eco:b3115 (2.7.2.15)
    eco:b2297 (2.3.1.8) + eco:b2463 (1.1.1.40)
    eco:b2320 (1.1.1.290) + eco:b1380 (1.1.1.28)
    eco:b2323 (2.3.1.41) + eco:b1095 (2.3.1.179)
    eco:b2373 (4.1.1.8) + eco:b0507 (4.1.1.47)
    eco:b2417 (2.7.1.-) + eco:b0679 (2.7.1.193)
    eco:b2429 (2.7.1.192) + eco:b2715 (2.7.1.-)
    eco:b2429 (2.7.1.192) + eco:b4240 (2.7.1.201)
    eco:b2463 (1.1.1.40) + eco:b2297 (2.3.1.8)
    eco:b2465 (2.2.1.1) + eco:b0420 (2.2.1.7)
    eco:b2478 (4.3.3.7) + eco:b3225 (4.1.3.3)
    eco:b2500 (2.1.2.2) + eco:b1232 (3.5.1.10)
    eco:b2507 (6.3.5.2) + eco:b3360 (2.6.1.85)
    eco:b2508 (1.1.1.205) + eco:b0104 (1.7.1.7)
    eco:b2523 (3.4.11.23) + eco:b4260 (3.4.11.1)
    eco:b2533 (3.1.3.25) + eco:b4214 (3.1.3.7)
    eco:b2541 (1.3.1.87) + eco:b2842 (1.1.1.127)
    eco:b2542 (1.18.1.3) + eco:b3365 (1.7.1.15)
    eco:b2574 (1.4.3.16) + eco:b4154 (1.3.5.4)
    eco:b2705 (1.1.1.140) + eco:b2842 (1.1.1.127)
    eco:b2715 (2.7.1.-) + eco:b2429 (2.7.1.192)
    eco:b2715 (2.7.1.-) + eco:b3599 (2.7.1.197)
    eco:b2715 (2.7.1.-) + eco:b4240 (2.7.1.201)
    eco:b2799 (1.1.1.77) + eco:b3011 (1.1.-.-)
    eco:b2800 (4.1.2.17) + eco:b4198 (5.1.3.4)
    eco:b2803 (2.7.1.51) + eco:b3564 (2.7.1.17)
    eco:b2803 (2.7.1.51) + eco:b3580 (2.7.1.53)
    eco:b2803 (2.7.1.51) + eco:b3904 (2.7.1.5)
    eco:b2842 (1.1.1.127) + eco:b0596 (1.3.1.28)
    eco:b2842 (1.1.1.127) + eco:b2541 (1.3.1.87)
    eco:b2842 (1.1.1.127) + eco:b2705 (1.1.1.140)
    eco:b2873 (3.5.2.2) + eco:b0512 (3.5.2.5)
    eco:b2878 (1.97.1.9) + eco:b2146 (1.3.1.1)
    eco:b2919 (4.1.1.41) + eco:b1393 (4.2.1.17)
    eco:b2919 (4.1.1.41) + eco:b1394 (5.3.3.18)
    eco:b2919 (4.1.1.41) + eco:b2262 (4.1.3.36)
    eco:b2925 (4.1.2.13) + eco:b2096 (4.1.2.40)
    eco:b2925 (4.1.2.13) + eco:b3137 (4.1.2.40)
    eco:b2927 (1.2.1.72) + eco:b1779 (1.2.1.12)
    eco:b2934 (2.7.1.197) + eco:b4195 (2.7.1.194)
    eco:b2934 (2.7.1.197) + eco:b4302 (2.7.1.200)
    eco:b2935 (2.2.1.1) + eco:b0420 (2.2.1.7)
    eco:b2965 (4.1.1.17) + eco:b0186 (4.1.1.18)
    eco:b2965 (4.1.1.17) + eco:b4117 (4.1.1.19)
    eco:b2965 (4.1.1.17) + eco:b4131 (4.1.1.18)
    eco:b3011 (1.1.-.-) + eco:b2799 (1.1.1.77)
    eco:b3073 (2.6.1.82) + eco:b0774 (2.6.1.62)
    eco:b3073 (2.6.1.82) + eco:b1748 (2.6.1.81)
    eco:b3091 (4.2.1.7) + eco:b3128 (4.2.1.42)
    eco:b3115 (2.7.2.15) + eco:b2296 (2.7.2.1)
    eco:b3126 (4.1.2.20) + eco:b2245 (4.1.2.53)
    eco:b3128 (4.2.1.42) + eco:b3091 (4.2.1.7)
    eco:b3133 (2.7.1.-) + eco:b1817 (2.7.1.191)
    eco:b3138 (2.7.1.-) + eco:b1817 (2.7.1.191)
    eco:b3141 (5.3.1.-) + eco:b0678 (3.5.99.6)
    eco:b3176 (5.4.2.10) + eco:b0688 (5.4.2.2)
    eco:b3176 (5.4.2.10) + eco:b2048 (5.4.2.8)
    eco:b3222 (2.7.1.60) + eco:b0394 (2.7.1.4)
    eco:b3222 (2.7.1.60) + eco:b1119 (2.7.1.59)
    eco:b3281 (1.1.1.25) + eco:b1692 (1.1.1.282)
    eco:b3365 (1.7.1.15) + eco:b2542 (1.18.1.3)
    eco:b3365 (1.7.1.15) + eco:b3500 (1.8.1.7)
    eco:b3365 (1.7.1.15) + eco:b3962 (1.6.1.1)
    eco:b3385 (3.1.3.18) + eco:b1317 (5.4.2.6)
    eco:b3386 (5.1.3.1) + eco:b4085 (5.1.3.-)
    eco:b3431 (3.2.1.196) + eco:b3432 (2.4.1.18)
    eco:b3432 (2.4.1.18) + eco:b3431 (3.2.1.196)
    eco:b3500 (1.8.1.7) + eco:b3365 (1.7.1.15)
    eco:b3500 (1.8.1.7) + eco:b3962 (1.6.1.1)
    eco:b3551 (1.-.-.-) + eco:b0894 (1.8.5.3)
    eco:b3551 (1.-.-.-) + eco:b1587 (1.97.1.9)
    eco:b3551 (1.-.-.-) + eco:b1588 (1.97.1.9)
    eco:b3551 (1.-.-.-) + eco:b2206 (1.7.99.-)
    eco:b3564 (2.7.1.17) + eco:b0063 (2.7.1.16)
    eco:b3564 (2.7.1.17) + eco:b2803 (2.7.1.51)
    eco:b3564 (2.7.1.17) + eco:b3580 (2.7.1.53)
    eco:b3564 (2.7.1.17) + eco:b3904 (2.7.1.5)
    eco:b3571 (3.2.1.1) + eco:b4239 (3.2.1.93)
    eco:b3575 (1.1.1.130) + eco:b0517 (1.1.1.350)
    eco:b3580 (2.7.1.53) + eco:b2803 (2.7.1.51)
    eco:b3580 (2.7.1.53) + eco:b3564 (2.7.1.17)
    eco:b3580 (2.7.1.53) + eco:b3904 (2.7.1.5)
    eco:b3588 (1.2.1.-) + eco:b0312 (1.2.1.8)
    eco:b3588 (1.2.1.-) + eco:b1300 (1.2.1.99)
    eco:b3588 (1.2.1.-) + eco:b1385 (1.2.1.39)
    eco:b3588 (1.2.1.-) + eco:b1444 (1.2.1.19)
    eco:b3588 (1.2.1.-) + eco:b1746 (1.2.1.71)
    eco:b3589 (1.1.1.1) + eco:b2799 (1.1.1.77)
    eco:b3589 (1.1.1.1) + eco:b3011 (1.1.-.-)
    eco:b3599 (2.7.1.197) + eco:b2169 (2.7.1.202)
    eco:b3599 (2.7.1.197) + eco:b2715 (2.7.1.-)
    eco:b3599 (2.7.1.197) + eco:b4195 (2.7.1.194)
    eco:b3600 (1.1.1.17) + eco:b4323 (1.1.1.57)
    eco:b3616 (1.1.1.103) + eco:b1580 (1.1.1.380)
    eco:b3616 (1.1.1.103) + eco:b2091 (1.1.1.251)
    eco:b3617 (2.3.1.29) + eco:b0776 (2.3.1.47)
    eco:b3671 (2.2.1.6) + eco:b0507 (4.1.1.47)
    eco:b3671 (2.2.1.6) + eco:b2264 (2.2.1.9)
    eco:b3671 (2.2.1.6) + eco:b2373 (4.1.1.8)
    eco:b3752 (2.7.1.15) + eco:b1723 (2.7.1.11)
    eco:b3752 (2.7.1.15) + eco:b2166 (2.7.1.83)
    eco:b3752 (2.7.1.15) + eco:b2168 (2.7.1.56)
    eco:b3771 (4.2.1.9) + eco:b1851 (4.2.1.12)
    eco:b3788 (4.2.1.46) + eco:b2052 (1.1.1.271)
    eco:b3788 (4.2.1.46) + eco:b2053 (4.2.1.47)
    eco:b3831 (2.4.2.3) + eco:b4384 (2.4.2.1)
    eco:b3870 (6.3.1.2) + eco:b1297 (6.3.1.11)
    eco:b3894 (1.2.1.2) + eco:b0894 (1.8.5.3)
    eco:b3894 (1.2.1.2) + eco:b1587 (1.97.1.9)
    eco:b3894 (1.2.1.2) + eco:b1588 (1.97.1.9)
    eco:b3894 (1.2.1.2) + eco:b2206 (1.7.99.-)
    eco:b3904 (2.7.1.5) + eco:b2803 (2.7.1.51)
    eco:b3904 (2.7.1.5) + eco:b3564 (2.7.1.17)
    eco:b3904 (2.7.1.5) + eco:b3580 (2.7.1.53)
    eco:b3926 (2.7.1.30) + eco:b2803 (2.7.1.51)
    eco:b3926 (2.7.1.30) + eco:b3564 (2.7.1.17)
    eco:b3926 (2.7.1.30) + eco:b3580 (2.7.1.53)
    eco:b3939 (2.5.1.48) + eco:b3008 (4.4.1.8)
    eco:b3960 (4.3.2.1) + eco:b4139 (4.3.1.1)
    eco:b3962 (1.6.1.1) + eco:b3365 (1.7.1.15)
    eco:b3962 (1.6.1.1) + eco:b3500 (1.8.1.7)
    eco:b4019 (2.1.1.13) + eco:b0261 (2.1.1.10)
    eco:b4069 (6.2.1.1) + eco:b0335 (6.2.1.17)
    eco:b4069 (6.2.1.1) + eco:b0586 (6.3.2.14)
    eco:b4084 (2.7.1.55) + eco:b0394 (2.7.1.4)
    eco:b4084 (2.7.1.55) + eco:b1119 (2.7.1.59)
    eco:b4085 (5.1.3.-) + eco:b3386 (5.1.3.1)
    eco:b4117 (4.1.1.19) + eco:b0186 (4.1.1.18)
    eco:b4117 (4.1.1.19) + eco:b0693 (4.1.1.17)
    eco:b4117 (4.1.1.19) + eco:b2965 (4.1.1.17)
    eco:b4117 (4.1.1.19) + eco:b4131 (4.1.1.18)
    eco:b4119 (3.2.1.22) + eco:b1734 (3.2.1.86)
    eco:b4122 (4.2.1.2) + eco:b3061 (4.2.1.32)
    eco:b4122 (4.2.1.2) + eco:b3062 (4.2.1.32)
    eco:b4131 (4.1.1.18) + eco:b0693 (4.1.1.17)
    eco:b4131 (4.1.1.18) + eco:b2965 (4.1.1.17)
    eco:b4131 (4.1.1.18) + eco:b4117 (4.1.1.19)
    eco:b4195 (2.7.1.194) + eco:b2934 (2.7.1.197)
    eco:b4195 (2.7.1.194) + eco:b3599 (2.7.1.197)
    eco:b4195 (2.7.1.194) + eco:b4302 (2.7.1.200)
    eco:b4198 (5.1.3.4) + eco:b2800 (4.1.2.17)
    eco:b4239 (3.2.1.93) + eco:b1309 (2.4.1.7)
    eco:b4240 (2.7.1.201) + eco:b2429 (2.7.1.192)
    eco:b4240 (2.7.1.201) + eco:b2715 (2.7.1.-)
    eco:b4260 (3.4.11.1) + eco:b2523 (3.4.11.23)
    eco:b4298 (4.3.3.7) + eco:b3225 (4.1.3.3)
    eco:b4302 (2.7.1.200) + eco:b2934 (2.7.1.197)
    eco:b4302 (2.7.1.200) + eco:b4195 (2.7.1.194)
    eco:b4323 (1.1.1.57) + eco:b1521 (1.1.1.58)
    eco:b4323 (1.1.1.57) + eco:b3600 (1.1.1.17)
    eco:b4384 (2.4.2.1) + eco:b3831 (2.4.2.3)
    eco:b4395 (5.4.2.12) + eco:b0638 (3.1.3.73)
    eco:b4478 (4.2.1.6) + eco:b1581 (4.2.1.8)
    eco:b4478 (4.2.1.6) + eco:b2247 (4.2.1.90)
    
    148/413 -> 35.8 percent neofunctionalisation of all new functions.

Conclusion
----------
Many more neofunctionalisations are reported when no gene duplication 'triangle' has to be proven: 35% instead of 10% of all new EC numbers. In reverse, this means that for about 20% of new functions, which do have paralogs, no ortholog in sso could be found.
This might be due to:
a) gene duplications that happened before the common ancestor parted into sso and eco, and have since been lost/mutated in sso. Keep in mind that, in this experiment, the common ancestor ist surrogated by sso, which is a questionable assumption.
b) gene duplications that happened after the common ancestor parted into sso and eco, and only happened in eco. This seems to be the better explanation, because there have presumably been many evolutionary events between the last common ancestor and eco.
Therefore, this method is not suitable for comparison of two individual species, because the 'older one' has to surrogate the common ancestor, because reconstructing a common ancestor's genes is a non-trivial task.
Shifting the scope from individual species, surrogating one for the common ancestor, to a group of species, comparing it to a super-group, might solve the problem.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph, SubstanceEnzymeGraph
from FEV_KEGG.KEGG import Database
import FEV_KEGG.KEGG.Organism


if __name__ == '__main__':
    
    #- Get all known pathways of eco and sso from KEGG.
    eco = FEV_KEGG.KEGG.Organism.Organism('eco')
    sso = FEV_KEGG.KEGG.Organism.Organism('sso')
    
    ecoPathways = eco.getMetabolicPathways()
    ssoPathways = sso.getMetabolicPathways()
    
    #- Convert each pathway into a substance-ecNumber graph. This is the incomplete metabolic network.
    ecoReactionGraph = SubstanceReactionGraph.fromPathway(ecoPathways)
    ssoReactionGraph = SubstanceReactionGraph.fromPathway(ssoPathways)
    
    ecoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ecoReactionGraph)
    ssoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ssoReactionGraph)
    
    ecoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ecoGeneGraph)
    ssoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ssoGeneGraph)
    
    #- Remove multifunctional enzymes, meaning enzymes associated with more than one EC number. Helps to reduce false gene duplications.
    ecoEnzymeGraph.removeMultifunctionalEnzymes()
    ssoEnzymeGraph.removeMultifunctionalEnzymes()
    
    ecoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ecoEnzymeGraph)
    ssoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ssoEnzymeGraph)
    
    #- Combine both species' networks to a consensus network, by INTERSECTION operation.
    intersectionEcGraph = ecoEcGraph.intersection(ssoEcGraph)
    
    #- Create the networks of EC numbers that only occur in eco.
    onlyEcoEcGraph = ecoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    
    
    
    output = []
    
    #- For each eco-only EC number, gather any eco gene which encodes an enzyme associated with this EC number.
    new_ec_numbers = onlyEcoEcGraph.getECs()
    neofunctionalised_ec_numbers = set()
    for ecNumber in new_ec_numbers:
        geneIDs = ecoEnzymeGraph.getGeneIDsForEcNumber(ecNumber)
        
        #- For each such gene, find homologs within eco (paralogs).
        for geneID in geneIDs:
            paralogs = Database.getParalogsOnlyGeneID(geneID)
                        
            #- If there exists at least one paralog with a different function, report a neofunctionalisation.
            paralogs_different_function = []
            
            for paralog in paralogs:
                paralog_enzyme = ecoEnzymeGraph.getEnzymeForGeneID(paralog)
                
                if paralog_enzyme is None: # can happen, when the paralogous gene encodes a non-enzymatic protein
                    continue
                elif len(paralog_enzyme.ecNumbers) > 0 and ecNumber not in paralog_enzyme.ecNumbers:
                    output.append(paralog.__str__() + " (" + ", ".join([x.__str__() for x in paralog_enzyme.ecNumbers]) + ") + " + geneID.__str__() + " (" + ecNumber.__str__() + ")")
                    
                    #- Report percentage neofunctionalisations of all new functions.
                    neofunctionalised_ec_numbers.add(ecNumber)                    
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
        
    print("\n" + str(len(neofunctionalised_ec_numbers)) + "/" + str(len(new_ec_numbers)) + " -> %2.1f percent neofunctionalisation of all new functions." % (len(neofunctionalised_ec_numbers)/len(new_ec_numbers)*100))
        