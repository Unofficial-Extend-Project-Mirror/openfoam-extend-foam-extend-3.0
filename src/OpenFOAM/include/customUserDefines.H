/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Some components of OpenFOAM have customizable compile configurations.  They
    look to this file for compiler defines.  If you wish to use one of these
    non-standard configurations, make changes below.

\*---------------------------------------------------------------------------*/

/* VectorN library defines - this determines the available sizes of VectorN
 * (and associated) objects, e.g. Vector8, Tensor5, etc.
 * The range must be within 2 and 255, including 2 and 255.  Uncomment any
 * combination of defines below.
 */

#define VECTORN_2 /*default*/
#define VECTORN_4 /*default*/
#define VECTORN_5 /*default*/
#define VECTORN_6 /*default*/
#define VECTORN_7 /*default*/
#define VECTORN_8 /*default*/
//#define VECTORN_9
//#define VECTORN_10
//#define VECTORN_11
//#define VECTORN_12
//#define VECTORN_13
//#define VECTORN_14
//#define VECTORN_15
//#define VECTORN_16
//#define VECTORN_17
//#define VECTORN_18
//#define VECTORN_19
//#define VECTORN_20
//#define VECTORN_21
//#define VECTORN_22
//#define VECTORN_23
//#define VECTORN_24
//#define VECTORN_25
//#define VECTORN_26
//#define VECTORN_27
//#define VECTORN_28
//#define VECTORN_29
//#define VECTORN_30
//#define VECTORN_31
//#define VECTORN_32
//#define VECTORN_33
//#define VECTORN_34
//#define VECTORN_35
//#define VECTORN_36
//#define VECTORN_37
//#define VECTORN_38
//#define VECTORN_39
//#define VECTORN_40
//#define VECTORN_41
//#define VECTORN_42
//#define VECTORN_43
//#define VECTORN_44
//#define VECTORN_45
//#define VECTORN_46
//#define VECTORN_47
//#define VECTORN_48
//#define VECTORN_49
//#define VECTORN_50
//#define VECTORN_51
//#define VECTORN_52
//#define VECTORN_53
//#define VECTORN_54
//#define VECTORN_55
//#define VECTORN_56
//#define VECTORN_57
//#define VECTORN_58
//#define VECTORN_59
//#define VECTORN_60
//#define VECTORN_61
//#define VECTORN_62
//#define VECTORN_63
//#define VECTORN_64
//#define VECTORN_65
//#define VECTORN_66
//#define VECTORN_67
//#define VECTORN_68
//#define VECTORN_69
//#define VECTORN_70
//#define VECTORN_71
//#define VECTORN_72
//#define VECTORN_73
//#define VECTORN_74
//#define VECTORN_75
//#define VECTORN_76
//#define VECTORN_77
//#define VECTORN_78
//#define VECTORN_79
//#define VECTORN_80
//#define VECTORN_81
//#define VECTORN_82
//#define VECTORN_83
//#define VECTORN_84
//#define VECTORN_85
//#define VECTORN_86
//#define VECTORN_87
//#define VECTORN_88
//#define VECTORN_89
//#define VECTORN_90
//#define VECTORN_91
//#define VECTORN_92
//#define VECTORN_93
//#define VECTORN_94
//#define VECTORN_95
//#define VECTORN_96
//#define VECTORN_97
//#define VECTORN_98
//#define VECTORN_99
//#define VECTORN_100
//#define VECTORN_101
//#define VECTORN_102
//#define VECTORN_103
//#define VECTORN_104
//#define VECTORN_105
//#define VECTORN_106
//#define VECTORN_107
//#define VECTORN_108
//#define VECTORN_109
//#define VECTORN_110
//#define VECTORN_111
//#define VECTORN_112
//#define VECTORN_113
//#define VECTORN_114
//#define VECTORN_115
//#define VECTORN_116
//#define VECTORN_117
//#define VECTORN_118
//#define VECTORN_119
//#define VECTORN_120
//#define VECTORN_121
//#define VECTORN_122
//#define VECTORN_123
//#define VECTORN_124
//#define VECTORN_125
//#define VECTORN_126
//#define VECTORN_127
//#define VECTORN_128
//#define VECTORN_129
//#define VECTORN_130
//#define VECTORN_131
//#define VECTORN_132
//#define VECTORN_133
//#define VECTORN_134
//#define VECTORN_135
//#define VECTORN_136
//#define VECTORN_137
//#define VECTORN_138
//#define VECTORN_139
//#define VECTORN_140
//#define VECTORN_141
//#define VECTORN_142
//#define VECTORN_143
//#define VECTORN_144
//#define VECTORN_145
//#define VECTORN_146
//#define VECTORN_147
//#define VECTORN_148
//#define VECTORN_149
//#define VECTORN_150
//#define VECTORN_151
//#define VECTORN_152
//#define VECTORN_153
//#define VECTORN_154
//#define VECTORN_155
//#define VECTORN_156
//#define VECTORN_157
//#define VECTORN_158
//#define VECTORN_159
//#define VECTORN_160
//#define VECTORN_161
//#define VECTORN_162
//#define VECTORN_163
//#define VECTORN_164
//#define VECTORN_165
//#define VECTORN_166
//#define VECTORN_167
//#define VECTORN_168
//#define VECTORN_169
//#define VECTORN_170
//#define VECTORN_171
//#define VECTORN_172
//#define VECTORN_173
//#define VECTORN_174
//#define VECTORN_175
//#define VECTORN_176
//#define VECTORN_177
//#define VECTORN_178
//#define VECTORN_179
//#define VECTORN_180
//#define VECTORN_181
//#define VECTORN_182
//#define VECTORN_183
//#define VECTORN_184
//#define VECTORN_185
//#define VECTORN_186
//#define VECTORN_187
//#define VECTORN_188
//#define VECTORN_189
//#define VECTORN_190
//#define VECTORN_191
//#define VECTORN_192
//#define VECTORN_193
//#define VECTORN_194
//#define VECTORN_195
//#define VECTORN_196
//#define VECTORN_197
//#define VECTORN_198
//#define VECTORN_199
//#define VECTORN_200
//#define VECTORN_201
//#define VECTORN_202
//#define VECTORN_203
//#define VECTORN_204
//#define VECTORN_205
//#define VECTORN_206
//#define VECTORN_207
//#define VECTORN_208
//#define VECTORN_209
//#define VECTORN_210
//#define VECTORN_211
//#define VECTORN_212
//#define VECTORN_213
//#define VECTORN_214
//#define VECTORN_215
//#define VECTORN_216
//#define VECTORN_217
//#define VECTORN_218
//#define VECTORN_219
//#define VECTORN_220
//#define VECTORN_221
//#define VECTORN_222
//#define VECTORN_223
//#define VECTORN_224
//#define VECTORN_225
//#define VECTORN_226
//#define VECTORN_227
//#define VECTORN_228
//#define VECTORN_229
//#define VECTORN_230
//#define VECTORN_231
//#define VECTORN_232
//#define VECTORN_233
//#define VECTORN_234
//#define VECTORN_235
//#define VECTORN_236
//#define VECTORN_237
//#define VECTORN_238
//#define VECTORN_239
//#define VECTORN_240
//#define VECTORN_241
//#define VECTORN_242
//#define VECTORN_243
//#define VECTORN_244
//#define VECTORN_245
//#define VECTORN_246
//#define VECTORN_247
//#define VECTORN_248
//#define VECTORN_249
//#define VECTORN_250
//#define VECTORN_251
//#define VECTORN_252
//#define VECTORN_253
//#define VECTORN_254
//#define VECTORN_255


/* VectorSpace compile configuration - this affects the degree to which
 * template metaprogramming is used in the definition of VectorSpace.
 * If you are using high VECTORN groupings, this option has a significant
 * effect on the performance of your compile, and your runtime execution.
 *
 *    Subject       |   Less metaprogramming    |   More metaprogramming
 *   ---------------+---------------------------+--------------------------
 *   Execution time |          Slower           |         Faster
 *   Compile time   |          Faster           |         Slower
 *   Library size   |          Smaller          |         Larger
 *
 *   Options are:
 *      VECTORSPACE_FORLOOPS - Minimal metaprogramming
 *      VECTORSPACE_PARTIAL  - Moderate metaprogramming *note - see below
 *      VECTORSPACE_FULL     - Full metaprogramming
 *
 *   * note: If using VECTORSPACE_PARTIAL, you must also define
 *      VECTORSPACE_UNROLLN as an integer between 2 and 255.  A lower number
 *      means more metaprogramming, and a higher number means less.  A value of
 *      255 is essentially the same as VECTORSPACE_FULL, whereas a value of
 *      1 is very close to VECTORSPACE_FORLOOPS. VECTORSPACE_UNROLLN should
 *      be less than the maximum instantiation depth for template classes to
 *      avoid compiler errors
 */

// UNCOMMENT ONLY ONE OF THESE

// #define VECTORSPACE_FORLOOPS
#define VECTORSPACE_PARTIAL /*default*/
// #define VECTORSPACE_FULL

/*
 * If VECTORSPACE_PARTIAL is defined, use a number between 2 and 255 below
 */
 
#define VECTORSPACE_UNROLLN 16 /*default is 16*/

#include "customUserDefineConditionals.H"
