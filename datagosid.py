# Great Snub Icosidodecahedron

import numpy as np

C0 = 0.0812642829572718262058790219454
C1 = 0.113690409611386992694315866713
C2 = 0.125868396751318911853846136829
C3 = 0.239558806362705904548162003542
C4 = 0.257356768647573545542485428757
C5 = 0.284923627010362671688163165899
C6 = 0.335147715949345479170923435882
C7 = 0.371047178258960538236801295470
C8 = 0.468878573956484550146522281604
C9 = 0.519102662895467357629282551587
C10 = 0.574706522312051383719085439424
C11 = 0.644971059646786269483128688416
C12 = 0.655970805269323209924964461370
C13 = 0.726235342604058095689007710362
C14 = 0.804026289905830029317445717487

# C0  = square-root of a root of the polynomial:
#     4096*(x^6] + 3072*(x^5] - 3584*(x^4] - 2048*(x^3] + 1312*(x^2] - 160*x + 1
# C1  = square-root of a root of the polynomial:
#     4096*(x^6] - 13312*(x^5] + 9216*(x^4] - 9472*(x^3] + 1872*(x^2] - 100*x + 1
# C2  = square-root of a root of the polynomial:
#     4096*(x^6] - 8192*(x^5] + 1792*(x^4] - 7488*(x^3] + 3456*(x^2] - 116*x + 1
# C3  = square-root of a root of the polynomial:
#     4096*(x^6] - 15360*(x^5] + 18944*(x^4] - 7168*(x^3] + 1024*(x^2] - 56*x + 1
# C4  = square-root of a root of the polynomial:
#     4096*(x^6] - 19456*(x^5] + 14592*(x^4] - 4736*(x^3] + 752*(x^2] - 48*x + 1
# C5  = square-root of a root of the polynomial:
#     4096*(x^6] - 12288*(x^5] - 768*(x^4] + 384*(x^3] + 272*(x^2] - 36*x + 1
# C6  = square-root of a root of the polynomial:
#     4096*(x^6] + 6144*(x^5] + 4352*(x^4] - 3456*(x^3] + 672*(x^2] - 48*x + 1
# C7  = square-root of a root of the polynomial:
#     4096*(x^6] - 12288*(x^5] + 15872*(x^4] - 6016*(x^3] + 912*(x^2] - 56*x + 1
# C8  = square-root of a root of the polynomial:
#     4096*(x^6] - 1024*(x^5] + 4096*(x^4] - 4672*(x^3] + 1392*(x^2] - 128*x + 1
# C9  = square-root of a root of the polynomial:
#     4096*(x^6] - 5120*(x^5] + 9472*(x^4] - 5888*(x^3] + 1216*(x^2] - 68*x + 1
# C10 = square-root of a root of the polynomial:
#     4096*(x^6] - 3072*(x^5] + 9728*(x^4] - 8960*(x^3] + 2944*(x^2] - 328*x + 1
# C11 = square-root of a root of the polynomial:
#     4096*(x^6] - 21504*(x^5] + 16384*(x^4] - 4672*(x^3] + 624*(x^2] - 40*x + 1
# C12 = square-root of a root of the polynomial:  4096*(x^6] - 18432*(x^5]
#     + 29184*(x^4] - 20160*(x^3] + 5728*(x^2] - 488*x + 1
# C13 = square-root of a root of the polynomial:
#     4096*(x^6] - 9728*(x^4] - 3072*(x^3] + 4256*(x^2] - 132*x + 1
# C14 = square-root of a root of the polynomial:  4096*(x^6] - 17408*(x^5]
#     + 28672*(x^4] - 21696*(x^3] + 6672*(x^2] - 416*x + 1

V0 = [-C7,   C6, -C11]
V1 = [C7,   C6,  C11]
V2 = [C7,  -C6, -C11]
V3 = [-C7,  -C6,  C11]
V4 = [C6, -C11,  -C7]
V5 = [-C6, -C11,   C7]
V6 = [-C6,  C11,  -C7]
V7 = [C6,  C11,   C7]
V8 = [-C11,  -C7,   C6]
V9 = [C11,  -C7,  -C6]
V10 = [C11,   C7,   C6]
V11 = [-C11,   C7,  -C6]
V12 = [C9,  C10,  -C4]
V13 = [-C9,  C10,   C4]
V14 = [-C9, -C10,  -C4]
V15 = [C9, -C10,   C4]
V16 = [-C10,  -C4,  -C9]
V17 = [C10,  -C4,   C9]
V18 = [C10,   C4,  -C9]
V19 = [-C10,   C4,   C9]
V20 = [C4,  -C9,  C10]
V21 = [-C4,  -C9, -C10]
V22 = [-C4,   C9,  C10]
V23 = [C4,   C9, -C10]
V24 = [-C13,   C3,   C5]
V25 = [C13,   C3,  -C5]
V26 = [C13,  -C3,   C5]
V27 = [-C13,  -C3,  -C5]
V28 = [C3,   C5, -C13]
V29 = [-C3,   C5,  C13]
V30 = [-C3,  -C5, -C13]
V31 = [C3,  -C5,  C13]
V32 = [C5, -C13,   C3]
V33 = [-C5, -C13,  -C3]
V34 = [-C5,  C13,   C3]
V35 = [C5,  C13,  -C3]
V36 = [C0,  C14,   C1]
V37 = [-C0,  C14,  -C1]
V38 = [-C0, -C14,   C1]
V39 = [C0, -C14,  -C1]
V40 = [C14,   C1,   C0]
V41 = [-C14,   C1,  -C0]
V42 = [-C14,  -C1,   C0]
V43 = [C14,  -C1,  -C0]
V44 = [C1,   C0,  C14]
V45 = [-C1,   C0, -C14]
V46 = [-C1,  -C0,  C14]
V47 = [C1,  -C0, -C14]
V48 = [-C8,  C12,  -C2]
V49 = [C8,  C12,   C2]
V50 = [C8, -C12,  -C2]
V51 = [-C8, -C12,   C2]
V52 = [-C12,  -C2,   C8]
V53 = [C12,  -C2,  -C8]
V54 = [C12,   C2,   C8]
V55 = [-C12,   C2,  -C8]
V56 = [C2,   C8,  C12]
V57 = [-C2,   C8, -C12]
V58 = [-C2,  -C8,  C12]
V59 = [C2,  -C8, -C12]

# Faces:
pentagram = np.zeros((12, 5, 3))
pentagram[0] = [V0, V36, V28, V48, V12]
pentagram[1] = [V1, V37, V29, V49, V13]
pentagram[2] = [V2, V38, V30, V50, V14]
pentagram[3] = [V3, V39, V31, V51, V15]
pentagram[4] = [V4, V40, V32, V53, V17]
pentagram[5] = [V5, V41, V33, V52, V16]
pentagram[6] = [V6, V42, V34, V55, V19]
pentagram[7] = [V7, V43, V35, V54, V18]
pentagram[8] = [V8, V44, V24, V58, V22]
pentagram[9] = [V9, V45, V25, V59, V23]
pentagram[10] = [V10, V46, V26, V56, V20]
pentagram[11] = [V11, V47, V27, V57, V21]

triangle = np.zeros((80, 3, 3))
triangle[0] = [V0, V2, V14]
triangle[1] = [V1, V3, V15]
triangle[2] = [V2, V0, V12]
triangle[3] = [V3, V1, V13]
triangle[4] = [V4, V5, V16]
triangle[5] = [V5, V4, V17]
triangle[6] = [V6, V7, V18]
triangle[7] = [V7, V6, V19]
triangle[8] = [V8, V11, V21]
triangle[9] = [V9, V10, V20]
triangle[10] = [V10, V9, V23]
triangle[11] = [V11, V8, V22]
triangle[12] = [V12, V48, V56]
triangle[13] = [V13, V49, V57]
triangle[14] = [V14, V50, V58]
triangle[15] = [V15, V51, V59]
triangle[16] = [V16, V52, V48]
triangle[17] = [V17, V53, V49]
triangle[18] = [V18, V54, V50]
triangle[19] = [V19, V55, V51]
triangle[20] = [V20, V56, V52]
triangle[21] = [V21, V57, V53]
triangle[22] = [V22, V58, V54]
triangle[23] = [V23, V59, V55]
triangle[24] = [V24, V44, V36]
triangle[25] = [V25, V45, V37]
triangle[26] = [V26, V46, V38]
triangle[27] = [V27, V47, V39]
triangle[28] = [V28, V36, V40]
triangle[29] = [V29, V37, V41]
triangle[30] = [V30, V38, V42]
triangle[31] = [V31, V39, V43]
triangle[32] = [V32, V40, V44]
triangle[33] = [V33, V41, V45]
triangle[34] = [V34, V42, V46]
triangle[35] = [V35, V43, V47]
triangle[36] = [V36, V0, V24]
triangle[37] = [V37, V1, V25]
triangle[38] = [V38, V2, V26]
triangle[39] = [V39, V3, V27]
triangle[40] = [V40, V4, V28]
triangle[41] = [V41, V5, V29]
triangle[42] = [V42, V6, V30]
triangle[43] = [V43, V7, V31]
triangle[44] = [V44, V8, V32]
triangle[45] = [V45, V9, V33]
triangle[46] = [V46, V10, V34]
triangle[47] = [V47, V11, V35]
triangle[48] = [V48, V28, V16]
triangle[49] = [V49, V29, V17]
triangle[50] = [V50, V30, V18]
triangle[51] = [V51, V31, V19]
triangle[52] = [V52, V33, V20]
triangle[53] = [V53, V32, V21]
triangle[54] = [V54, V35, V22]
triangle[55] = [V55, V34, V23]
triangle[56] = [V56, V26, V12]
triangle[57] = [V57, V27, V13]
triangle[58] = [V58, V24, V14]
triangle[59] = [V59, V25, V15]
triangle[60] = [V24, V0, V14]
triangle[61] = [V25, V1, V15]
triangle[62] = [V26, V2, V12]
triangle[63] = [V27, V3, V13]
triangle[64] = [V28, V4, V16]
triangle[65] = [V29, V5, V17]
triangle[66] = [V30, V6, V18]
triangle[67] = [V31, V7, V19]
triangle[68] = [V32, V8, V21]
triangle[69] = [V33, V9, V20]
triangle[70] = [V34, V10, V23]
triangle[71] = [V35, V11, V22]
triangle[72] = [V36, V44, V40]
triangle[73] = [V37, V45, V41]
triangle[74] = [V38, V46, V42]
triangle[75] = [V39, V47, V43]
triangle[76] = [V48, V52, V56]
triangle[77] = [V49, V53, V57]
triangle[78] = [V50, V54, V58]
triangle[79] = [V51, V55, V59]

face = [triangle, pentagram]
