#include "sparse_multiplier.h"

__device__
void sparse_multiplier(const Real * A, const Real * Vm, Real* w) {
  w[0] = A[0] * Vm[0] +  A[14] * Vm[1] +  A[28] * Vm[2] +  A[42] * Vm[3] +  A[56] * Vm[4] +  A[70] * Vm[5] +  A[84] * Vm[6] +  A[98] * Vm[7] +  A[112] * Vm[8] +  A[126] * Vm[9] +  A[140] * Vm[10] +  A[154] * Vm[11] +  A[168] * Vm[12] +  A[182] * Vm[13];
  w[1] = A[1] * Vm[0] +  A[15] * Vm[1] +  A[29] * Vm[2] +  A[43] * Vm[3] +  A[57] * Vm[4] +  A[71] * Vm[5] +  A[85] * Vm[6] +  A[99] * Vm[7] +  A[113] * Vm[8] +  A[127] * Vm[9] +  A[141] * Vm[10] +  A[155] * Vm[11] +  A[169] * Vm[12] +  A[183] * Vm[13];
  w[2] = A[2] * Vm[0] +  A[16] * Vm[1] +  A[30] * Vm[2] +  A[44] * Vm[3] +  A[58] * Vm[4] +  A[72] * Vm[5] +  A[86] * Vm[6] +  A[100] * Vm[7] +  A[114] * Vm[8] +  A[128] * Vm[9] +  A[142] * Vm[10] +  A[156] * Vm[11] +  A[170] * Vm[12] +  A[184] * Vm[13];
  w[3] = A[3] * Vm[0] +  A[17] * Vm[1] +  A[31] * Vm[2] +  A[45] * Vm[3] +  A[59] * Vm[4] +  A[73] * Vm[5] +  A[87] * Vm[6] +  A[101] * Vm[7] +  A[115] * Vm[8] +  A[129] * Vm[9] +  A[143] * Vm[10] +  A[157] * Vm[11] +  A[171] * Vm[12] +  A[185] * Vm[13];
  w[4] = A[4] * Vm[0] +  A[18] * Vm[1] +  A[32] * Vm[2] +  A[46] * Vm[3] +  A[60] * Vm[4] +  A[74] * Vm[5] +  A[88] * Vm[6] +  A[102] * Vm[7] +  A[116] * Vm[8] +  A[130] * Vm[9] +  A[144] * Vm[10] +  A[158] * Vm[11] +  A[172] * Vm[12] +  A[186] * Vm[13];
  w[5] = A[5] * Vm[0] +  A[19] * Vm[1] +  A[33] * Vm[2] +  A[47] * Vm[3] +  A[61] * Vm[4] +  A[75] * Vm[5] +  A[89] * Vm[6] +  A[103] * Vm[7] +  A[117] * Vm[8] +  A[131] * Vm[9] +  A[145] * Vm[10] +  A[159] * Vm[11] +  A[173] * Vm[12] +  A[187] * Vm[13];
  w[6] = A[6] * Vm[0] +  A[20] * Vm[1] +  A[34] * Vm[2] +  A[48] * Vm[3] +  A[62] * Vm[4] +  A[76] * Vm[5] +  A[90] * Vm[6] +  A[104] * Vm[7] +  A[118] * Vm[8] +  A[132] * Vm[9] +  A[146] * Vm[10] +  A[160] * Vm[11] +  A[174] * Vm[12] +  A[188] * Vm[13];
  w[7] = A[7] * Vm[0] +  A[21] * Vm[1] +  A[35] * Vm[2] +  A[49] * Vm[3] +  A[63] * Vm[4] +  A[77] * Vm[5] +  A[91] * Vm[6] +  A[105] * Vm[7] +  A[119] * Vm[8] +  A[133] * Vm[9] +  A[147] * Vm[10] +  A[161] * Vm[11] +  A[175] * Vm[12] +  A[189] * Vm[13];
  w[8] = A[8] * Vm[0] +  A[22] * Vm[1] +  A[36] * Vm[2] +  A[50] * Vm[3] +  A[64] * Vm[4] +  A[78] * Vm[5] +  A[92] * Vm[6] +  A[106] * Vm[7] +  A[120] * Vm[8] +  A[134] * Vm[9] +  A[148] * Vm[10] +  A[162] * Vm[11] +  A[176] * Vm[12] +  A[190] * Vm[13];
}
