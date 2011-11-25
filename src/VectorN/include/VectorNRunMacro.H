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
    A macro-runner for the VectorN library.  Executes macros only for active
    instantiations in the VectorN library.
    
    Usage:
    @verbatim
        
        #define VectorNMacro macroName
        #include "VectorNRunMacro.H"
        // VectorNMacro is automatically #undef'd
        
        // Or with optional arguments
        #define VectorNMacro macroName
        #define MacroArgs macroArguments
        #include "VectorNRunMacro.H"
        // VectorNMacro and MacroArgs are automatically undef'd
        
        // To omit a specific VECTORN_## (regardless of whether it exists or
        // not), use:
        #define SKIP_VECTORN_##
        // insert the lines that run the macros here (from above)
        #undef SKIP_VECTORN_##

        // This will have the effect of executing the macroName for only those
        // numbers that correspond to active VECTORN_## instantiations
        // e.g.:
        //
        // macroName(5);
        // macroName(8);
        // macroName(10);
        //
        // Or, when MacroArgs are included:
        //
        // macroName(5, macroArguments);
        // macroName(8, macroArguments);
        // macroName(10, macroArguments);
        
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef VectorNMacro
#   error VectorNMacro must be defined when including VectorNRunMacro.H
#endif

#ifdef MacroArgs
#   define RunMe(index) VectorNMacro(index, MacroArgs)
#else
#   define RunMe(index) VectorNMacro(index)
#endif

// There is no VECTORN_1 - use a scalar
#if defined (VECTORN_2) && !defined (SKIP_VECTORN_2)
    RunMe(2);
#endif
// There is no VECTORN_3 - use a vector
#if defined (VECTORN_4) && !defined (SKIP_VECTORN_4)
    RunMe(4);
#endif
#if defined (VECTORN_5) && !defined (SKIP_VECTORN_5)
    RunMe(5);
#endif
#if defined (VECTORN_6) && !defined (SKIP_VECTORN_6)
    RunMe(6);
#endif
#if defined (VECTORN_7) && !defined (SKIP_VECTORN_7)
    RunMe(7);
#endif
#if defined (VECTORN_8) && !defined (SKIP_VECTORN_8)
    RunMe(8);
#endif
#if defined (VECTORN_9) && !defined (SKIP_VECTORN_9)
    RunMe(9);
#endif
#if defined (VECTORN_10) && !defined (SKIP_VECTORN_10)
    RunMe(10);
#endif
#if defined (VECTORN_11) && !defined (SKIP_VECTORN_11)
    RunMe(11);
#endif
#if defined (VECTORN_12) && !defined (SKIP_VECTORN_12)
    RunMe(12);
#endif
#if defined (VECTORN_13) && !defined (SKIP_VECTORN_13)
    RunMe(13);
#endif
#if defined (VECTORN_14) && !defined (SKIP_VECTORN_14)
    RunMe(14);
#endif
#if defined (VECTORN_15) && !defined (SKIP_VECTORN_15)
    RunMe(15);
#endif
#if defined (VECTORN_16) && !defined (SKIP_VECTORN_16)
    RunMe(16);
#endif
#if defined (VECTORN_17) && !defined (SKIP_VECTORN_17)
    RunMe(17);
#endif
#if defined (VECTORN_18) && !defined (SKIP_VECTORN_18)
    RunMe(18);
#endif
#if defined (VECTORN_19) && !defined (SKIP_VECTORN_19)
    RunMe(19);
#endif
#if defined (VECTORN_20) && !defined (SKIP_VECTORN_20)
    RunMe(20);
#endif
#if defined (VECTORN_21) && !defined (SKIP_VECTORN_21)
    RunMe(21);
#endif
#if defined (VECTORN_22) && !defined (SKIP_VECTORN_22)
    RunMe(22);
#endif
#if defined (VECTORN_23) && !defined (SKIP_VECTORN_23)
    RunMe(23);
#endif
#if defined (VECTORN_24) && !defined (SKIP_VECTORN_24)
    RunMe(24);
#endif
#if defined (VECTORN_25) && !defined (SKIP_VECTORN_25)
    RunMe(25);
#endif
#if defined (VECTORN_26) && !defined (SKIP_VECTORN_26)
    RunMe(26);
#endif
#if defined (VECTORN_27) && !defined (SKIP_VECTORN_27)
    RunMe(27);
#endif
#if defined (VECTORN_28) && !defined (SKIP_VECTORN_28)
    RunMe(28);
#endif
#if defined (VECTORN_29) && !defined (SKIP_VECTORN_29)
    RunMe(29);
#endif
#if defined (VECTORN_30) && !defined (SKIP_VECTORN_30)
    RunMe(30);
#endif
#if defined (VECTORN_31) && !defined (SKIP_VECTORN_31)
    RunMe(31);
#endif
#if defined (VECTORN_32) && !defined (SKIP_VECTORN_32)
    RunMe(32);
#endif
#if defined (VECTORN_33) && !defined (SKIP_VECTORN_33)
    RunMe(33);
#endif
#if defined (VECTORN_34) && !defined (SKIP_VECTORN_34)
    RunMe(34);
#endif
#if defined (VECTORN_35) && !defined (SKIP_VECTORN_35)
    RunMe(35);
#endif
#if defined (VECTORN_36) && !defined (SKIP_VECTORN_36)
    RunMe(36);
#endif
#if defined (VECTORN_37) && !defined (SKIP_VECTORN_37)
    RunMe(37);
#endif
#if defined (VECTORN_38) && !defined (SKIP_VECTORN_38)
    RunMe(38);
#endif
#if defined (VECTORN_39) && !defined (SKIP_VECTORN_39)
    RunMe(39);
#endif
#if defined (VECTORN_40) && !defined (SKIP_VECTORN_40)
    RunMe(40);
#endif
#if defined (VECTORN_41) && !defined (SKIP_VECTORN_41)
    RunMe(41);
#endif
#if defined (VECTORN_42) && !defined (SKIP_VECTORN_42)
    RunMe(42);
#endif
#if defined (VECTORN_43) && !defined (SKIP_VECTORN_43)
    RunMe(43);
#endif
#if defined (VECTORN_44) && !defined (SKIP_VECTORN_44)
    RunMe(44);
#endif
#if defined (VECTORN_45) && !defined (SKIP_VECTORN_45)
    RunMe(45);
#endif
#if defined (VECTORN_46) && !defined (SKIP_VECTORN_46)
    RunMe(46);
#endif
#if defined (VECTORN_47) && !defined (SKIP_VECTORN_47)
    RunMe(47);
#endif
#if defined (VECTORN_48) && !defined (SKIP_VECTORN_48)
    RunMe(48);
#endif
#if defined (VECTORN_49) && !defined (SKIP_VECTORN_49)
    RunMe(49);
#endif
#if defined (VECTORN_50) && !defined (SKIP_VECTORN_50)
    RunMe(50);
#endif
#if defined (VECTORN_51) && !defined (SKIP_VECTORN_51)
    RunMe(51);
#endif
#if defined (VECTORN_52) && !defined (SKIP_VECTORN_52)
    RunMe(52);
#endif
#if defined (VECTORN_53) && !defined (SKIP_VECTORN_53)
    RunMe(53);
#endif
#if defined (VECTORN_54) && !defined (SKIP_VECTORN_54)
    RunMe(54);
#endif
#if defined (VECTORN_55) && !defined (SKIP_VECTORN_55)
    RunMe(55);
#endif
#if defined (VECTORN_56) && !defined (SKIP_VECTORN_56)
    RunMe(56);
#endif
#if defined (VECTORN_57) && !defined (SKIP_VECTORN_57)
    RunMe(57);
#endif
#if defined (VECTORN_58) && !defined (SKIP_VECTORN_58)
    RunMe(58);
#endif
#if defined (VECTORN_59) && !defined (SKIP_VECTORN_59)
    RunMe(59);
#endif
#if defined (VECTORN_60) && !defined (SKIP_VECTORN_60)
    RunMe(60);
#endif
#if defined (VECTORN_61) && !defined (SKIP_VECTORN_61)
    RunMe(61);
#endif
#if defined (VECTORN_62) && !defined (SKIP_VECTORN_62)
    RunMe(62);
#endif
#if defined (VECTORN_63) && !defined (SKIP_VECTORN_63)
    RunMe(63);
#endif
#if defined (VECTORN_64) && !defined (SKIP_VECTORN_64)
    RunMe(64);
#endif
#if defined (VECTORN_65) && !defined (SKIP_VECTORN_65)
    RunMe(65);
#endif
#if defined (VECTORN_66) && !defined (SKIP_VECTORN_66)
    RunMe(66);
#endif
#if defined (VECTORN_67) && !defined (SKIP_VECTORN_67)
    RunMe(67);
#endif
#if defined (VECTORN_68) && !defined (SKIP_VECTORN_68)
    RunMe(68);
#endif
#if defined (VECTORN_69) && !defined (SKIP_VECTORN_69)
    RunMe(69);
#endif
#if defined (VECTORN_70) && !defined (SKIP_VECTORN_70)
    RunMe(70);
#endif
#if defined (VECTORN_71) && !defined (SKIP_VECTORN_71)
    RunMe(71);
#endif
#if defined (VECTORN_72) && !defined (SKIP_VECTORN_72)
    RunMe(72);
#endif
#if defined (VECTORN_73) && !defined (SKIP_VECTORN_73)
    RunMe(73);
#endif
#if defined (VECTORN_74) && !defined (SKIP_VECTORN_74)
    RunMe(74);
#endif
#if defined (VECTORN_75) && !defined (SKIP_VECTORN_75)
    RunMe(75);
#endif
#if defined (VECTORN_76) && !defined (SKIP_VECTORN_76)
    RunMe(76);
#endif
#if defined (VECTORN_77) && !defined (SKIP_VECTORN_77)
    RunMe(77);
#endif
#if defined (VECTORN_78) && !defined (SKIP_VECTORN_78)
    RunMe(78);
#endif
#if defined (VECTORN_79) && !defined (SKIP_VECTORN_79)
    RunMe(79);
#endif
#if defined (VECTORN_80) && !defined (SKIP_VECTORN_80)
    RunMe(80);
#endif
#if defined (VECTORN_81) && !defined (SKIP_VECTORN_81)
    RunMe(81);
#endif
#if defined (VECTORN_82) && !defined (SKIP_VECTORN_82)
    RunMe(82);
#endif
#if defined (VECTORN_83) && !defined (SKIP_VECTORN_83)
    RunMe(83);
#endif
#if defined (VECTORN_84) && !defined (SKIP_VECTORN_84)
    RunMe(84);
#endif
#if defined (VECTORN_85) && !defined (SKIP_VECTORN_85)
    RunMe(85);
#endif
#if defined (VECTORN_86) && !defined (SKIP_VECTORN_86)
    RunMe(86);
#endif
#if defined (VECTORN_87) && !defined (SKIP_VECTORN_87)
    RunMe(87);
#endif
#if defined (VECTORN_88) && !defined (SKIP_VECTORN_88)
    RunMe(88);
#endif
#if defined (VECTORN_89) && !defined (SKIP_VECTORN_89)
    RunMe(89);
#endif
#if defined (VECTORN_90) && !defined (SKIP_VECTORN_90)
    RunMe(90);
#endif
#if defined (VECTORN_91) && !defined (SKIP_VECTORN_91)
    RunMe(91);
#endif
#if defined (VECTORN_92) && !defined (SKIP_VECTORN_92)
    RunMe(92);
#endif
#if defined (VECTORN_93) && !defined (SKIP_VECTORN_93)
    RunMe(93);
#endif
#if defined (VECTORN_94) && !defined (SKIP_VECTORN_94)
    RunMe(94);
#endif
#if defined (VECTORN_95) && !defined (SKIP_VECTORN_95)
    RunMe(95);
#endif
#if defined (VECTORN_96) && !defined (SKIP_VECTORN_96)
    RunMe(96);
#endif
#if defined (VECTORN_97) && !defined (SKIP_VECTORN_97)
    RunMe(97);
#endif
#if defined (VECTORN_98) && !defined (SKIP_VECTORN_98)
    RunMe(98);
#endif
#if defined (VECTORN_99) && !defined (SKIP_VECTORN_99)
    RunMe(99);
#endif
#if defined (VECTORN_100) && !defined (SKIP_VECTORN_100)
    RunMe(100);
#endif
#if defined (VECTORN_101) && !defined (SKIP_VECTORN_101)
    RunMe(101);
#endif
#if defined (VECTORN_102) && !defined (SKIP_VECTORN_102)
    RunMe(102);
#endif
#if defined (VECTORN_103) && !defined (SKIP_VECTORN_103)
    RunMe(103);
#endif
#if defined (VECTORN_104) && !defined (SKIP_VECTORN_104)
    RunMe(104);
#endif
#if defined (VECTORN_105) && !defined (SKIP_VECTORN_105)
    RunMe(105);
#endif
#if defined (VECTORN_106) && !defined (SKIP_VECTORN_106)
    RunMe(106);
#endif
#if defined (VECTORN_107) && !defined (SKIP_VECTORN_107)
    RunMe(107);
#endif
#if defined (VECTORN_108) && !defined (SKIP_VECTORN_108)
    RunMe(108);
#endif
#if defined (VECTORN_109) && !defined (SKIP_VECTORN_109)
    RunMe(109);
#endif
#if defined (VECTORN_110) && !defined (SKIP_VECTORN_110)
    RunMe(110);
#endif
#if defined (VECTORN_111) && !defined (SKIP_VECTORN_111)
    RunMe(111);
#endif
#if defined (VECTORN_112) && !defined (SKIP_VECTORN_112)
    RunMe(112);
#endif
#if defined (VECTORN_113) && !defined (SKIP_VECTORN_113)
    RunMe(113);
#endif
#if defined (VECTORN_114) && !defined (SKIP_VECTORN_114)
    RunMe(114);
#endif
#if defined (VECTORN_115) && !defined (SKIP_VECTORN_115)
    RunMe(115);
#endif
#if defined (VECTORN_116) && !defined (SKIP_VECTORN_116)
    RunMe(116);
#endif
#if defined (VECTORN_117) && !defined (SKIP_VECTORN_117)
    RunMe(117);
#endif
#if defined (VECTORN_118) && !defined (SKIP_VECTORN_118)
    RunMe(118);
#endif
#if defined (VECTORN_119) && !defined (SKIP_VECTORN_119)
    RunMe(119);
#endif
#if defined (VECTORN_120) && !defined (SKIP_VECTORN_120)
    RunMe(120);
#endif
#if defined (VECTORN_121) && !defined (SKIP_VECTORN_121)
    RunMe(121);
#endif
#if defined (VECTORN_122) && !defined (SKIP_VECTORN_122)
    RunMe(122);
#endif
#if defined (VECTORN_123) && !defined (SKIP_VECTORN_123)
    RunMe(123);
#endif
#if defined (VECTORN_124) && !defined (SKIP_VECTORN_124)
    RunMe(124);
#endif
#if defined (VECTORN_125) && !defined (SKIP_VECTORN_125)
    RunMe(125);
#endif
#if defined (VECTORN_126) && !defined (SKIP_VECTORN_126)
    RunMe(126);
#endif
#if defined (VECTORN_127) && !defined (SKIP_VECTORN_127)
    RunMe(127);
#endif
#if defined (VECTORN_128) && !defined (SKIP_VECTORN_128)
    RunMe(128);
#endif
#if defined (VECTORN_129) && !defined (SKIP_VECTORN_129)
    RunMe(129);
#endif
#if defined (VECTORN_130) && !defined (SKIP_VECTORN_130)
    RunMe(130);
#endif
#if defined (VECTORN_131) && !defined (SKIP_VECTORN_131)
    RunMe(131);
#endif
#if defined (VECTORN_132) && !defined (SKIP_VECTORN_132)
    RunMe(132);
#endif
#if defined (VECTORN_133) && !defined (SKIP_VECTORN_133)
    RunMe(133);
#endif
#if defined (VECTORN_134) && !defined (SKIP_VECTORN_134)
    RunMe(134);
#endif
#if defined (VECTORN_135) && !defined (SKIP_VECTORN_135)
    RunMe(135);
#endif
#if defined (VECTORN_136) && !defined (SKIP_VECTORN_136)
    RunMe(136);
#endif
#if defined (VECTORN_137) && !defined (SKIP_VECTORN_137)
    RunMe(137);
#endif
#if defined (VECTORN_138) && !defined (SKIP_VECTORN_138)
    RunMe(138);
#endif
#if defined (VECTORN_139) && !defined (SKIP_VECTORN_139)
    RunMe(139);
#endif
#if defined (VECTORN_140) && !defined (SKIP_VECTORN_140)
    RunMe(140);
#endif
#if defined (VECTORN_141) && !defined (SKIP_VECTORN_141)
    RunMe(141);
#endif
#if defined (VECTORN_142) && !defined (SKIP_VECTORN_142)
    RunMe(142);
#endif
#if defined (VECTORN_143) && !defined (SKIP_VECTORN_143)
    RunMe(143);
#endif
#if defined (VECTORN_144) && !defined (SKIP_VECTORN_144)
    RunMe(144);
#endif
#if defined (VECTORN_145) && !defined (SKIP_VECTORN_145)
    RunMe(145);
#endif
#if defined (VECTORN_146) && !defined (SKIP_VECTORN_146)
    RunMe(146);
#endif
#if defined (VECTORN_147) && !defined (SKIP_VECTORN_147)
    RunMe(147);
#endif
#if defined (VECTORN_148) && !defined (SKIP_VECTORN_148)
    RunMe(148);
#endif
#if defined (VECTORN_149) && !defined (SKIP_VECTORN_149)
    RunMe(149);
#endif
#if defined (VECTORN_150) && !defined (SKIP_VECTORN_150)
    RunMe(150);
#endif
#if defined (VECTORN_151) && !defined (SKIP_VECTORN_151)
    RunMe(151);
#endif
#if defined (VECTORN_152) && !defined (SKIP_VECTORN_152)
    RunMe(152);
#endif
#if defined (VECTORN_153) && !defined (SKIP_VECTORN_153)
    RunMe(153);
#endif
#if defined (VECTORN_154) && !defined (SKIP_VECTORN_154)
    RunMe(154);
#endif
#if defined (VECTORN_155) && !defined (SKIP_VECTORN_155)
    RunMe(155);
#endif
#if defined (VECTORN_156) && !defined (SKIP_VECTORN_156)
    RunMe(156);
#endif
#if defined (VECTORN_157) && !defined (SKIP_VECTORN_157)
    RunMe(157);
#endif
#if defined (VECTORN_158) && !defined (SKIP_VECTORN_158)
    RunMe(158);
#endif
#if defined (VECTORN_159) && !defined (SKIP_VECTORN_159)
    RunMe(159);
#endif
#if defined (VECTORN_160) && !defined (SKIP_VECTORN_160)
    RunMe(160);
#endif
#if defined (VECTORN_161) && !defined (SKIP_VECTORN_161)
    RunMe(161);
#endif
#if defined (VECTORN_162) && !defined (SKIP_VECTORN_162)
    RunMe(162);
#endif
#if defined (VECTORN_163) && !defined (SKIP_VECTORN_163)
    RunMe(163);
#endif
#if defined (VECTORN_164) && !defined (SKIP_VECTORN_164)
    RunMe(164);
#endif
#if defined (VECTORN_165) && !defined (SKIP_VECTORN_165)
    RunMe(165);
#endif
#if defined (VECTORN_166) && !defined (SKIP_VECTORN_166)
    RunMe(166);
#endif
#if defined (VECTORN_167) && !defined (SKIP_VECTORN_167)
    RunMe(167);
#endif
#if defined (VECTORN_168) && !defined (SKIP_VECTORN_168)
    RunMe(168);
#endif
#if defined (VECTORN_169) && !defined (SKIP_VECTORN_169)
    RunMe(169);
#endif
#if defined (VECTORN_170) && !defined (SKIP_VECTORN_170)
    RunMe(170);
#endif
#if defined (VECTORN_171) && !defined (SKIP_VECTORN_171)
    RunMe(171);
#endif
#if defined (VECTORN_172) && !defined (SKIP_VECTORN_172)
    RunMe(172);
#endif
#if defined (VECTORN_173) && !defined (SKIP_VECTORN_173)
    RunMe(173);
#endif
#if defined (VECTORN_174) && !defined (SKIP_VECTORN_174)
    RunMe(174);
#endif
#if defined (VECTORN_175) && !defined (SKIP_VECTORN_175)
    RunMe(175);
#endif
#if defined (VECTORN_176) && !defined (SKIP_VECTORN_176)
    RunMe(176);
#endif
#if defined (VECTORN_177) && !defined (SKIP_VECTORN_177)
    RunMe(177);
#endif
#if defined (VECTORN_178) && !defined (SKIP_VECTORN_178)
    RunMe(178);
#endif
#if defined (VECTORN_179) && !defined (SKIP_VECTORN_179)
    RunMe(179);
#endif
#if defined (VECTORN_180) && !defined (SKIP_VECTORN_180)
    RunMe(180);
#endif
#if defined (VECTORN_181) && !defined (SKIP_VECTORN_181)
    RunMe(181);
#endif
#if defined (VECTORN_182) && !defined (SKIP_VECTORN_182)
    RunMe(182);
#endif
#if defined (VECTORN_183) && !defined (SKIP_VECTORN_183)
    RunMe(183);
#endif
#if defined (VECTORN_184) && !defined (SKIP_VECTORN_184)
    RunMe(184);
#endif
#if defined (VECTORN_185) && !defined (SKIP_VECTORN_185)
    RunMe(185);
#endif
#if defined (VECTORN_186) && !defined (SKIP_VECTORN_186)
    RunMe(186);
#endif
#if defined (VECTORN_187) && !defined (SKIP_VECTORN_187)
    RunMe(187);
#endif
#if defined (VECTORN_188) && !defined (SKIP_VECTORN_188)
    RunMe(188);
#endif
#if defined (VECTORN_189) && !defined (SKIP_VECTORN_189)
    RunMe(189);
#endif
#if defined (VECTORN_190) && !defined (SKIP_VECTORN_190)
    RunMe(190);
#endif
#if defined (VECTORN_191) && !defined (SKIP_VECTORN_191)
    RunMe(191);
#endif
#if defined (VECTORN_192) && !defined (SKIP_VECTORN_192)
    RunMe(192);
#endif
#if defined (VECTORN_193) && !defined (SKIP_VECTORN_193)
    RunMe(193);
#endif
#if defined (VECTORN_194) && !defined (SKIP_VECTORN_194)
    RunMe(194);
#endif
#if defined (VECTORN_195) && !defined (SKIP_VECTORN_195)
    RunMe(195);
#endif
#if defined (VECTORN_196) && !defined (SKIP_VECTORN_196)
    RunMe(196);
#endif
#if defined (VECTORN_197) && !defined (SKIP_VECTORN_197)
    RunMe(197);
#endif
#if defined (VECTORN_198) && !defined (SKIP_VECTORN_198)
    RunMe(198);
#endif
#if defined (VECTORN_199) && !defined (SKIP_VECTORN_199)
    RunMe(199);
#endif
#if defined (VECTORN_200) && !defined (SKIP_VECTORN_200)
    RunMe(200);
#endif
#if defined (VECTORN_201) && !defined (SKIP_VECTORN_201)
    RunMe(201);
#endif
#if defined (VECTORN_202) && !defined (SKIP_VECTORN_202)
    RunMe(202);
#endif
#if defined (VECTORN_203) && !defined (SKIP_VECTORN_203)
    RunMe(203);
#endif
#if defined (VECTORN_204) && !defined (SKIP_VECTORN_204)
    RunMe(204);
#endif
#if defined (VECTORN_205) && !defined (SKIP_VECTORN_205)
    RunMe(205);
#endif
#if defined (VECTORN_206) && !defined (SKIP_VECTORN_206)
    RunMe(206);
#endif
#if defined (VECTORN_207) && !defined (SKIP_VECTORN_207)
    RunMe(207);
#endif
#if defined (VECTORN_208) && !defined (SKIP_VECTORN_208)
    RunMe(208);
#endif
#if defined (VECTORN_209) && !defined (SKIP_VECTORN_209)
    RunMe(209);
#endif
#if defined (VECTORN_210) && !defined (SKIP_VECTORN_210)
    RunMe(210);
#endif
#if defined (VECTORN_211) && !defined (SKIP_VECTORN_211)
    RunMe(211);
#endif
#if defined (VECTORN_212) && !defined (SKIP_VECTORN_212)
    RunMe(212);
#endif
#if defined (VECTORN_213) && !defined (SKIP_VECTORN_213)
    RunMe(213);
#endif
#if defined (VECTORN_214) && !defined (SKIP_VECTORN_214)
    RunMe(214);
#endif
#if defined (VECTORN_215) && !defined (SKIP_VECTORN_215)
    RunMe(215);
#endif
#if defined (VECTORN_216) && !defined (SKIP_VECTORN_216)
    RunMe(216);
#endif
#if defined (VECTORN_217) && !defined (SKIP_VECTORN_217)
    RunMe(217);
#endif
#if defined (VECTORN_218) && !defined (SKIP_VECTORN_218)
    RunMe(218);
#endif
#if defined (VECTORN_219) && !defined (SKIP_VECTORN_219)
    RunMe(219);
#endif
#if defined (VECTORN_220) && !defined (SKIP_VECTORN_220)
    RunMe(220);
#endif
#if defined (VECTORN_221) && !defined (SKIP_VECTORN_221)
    RunMe(221);
#endif
#if defined (VECTORN_222) && !defined (SKIP_VECTORN_222)
    RunMe(222);
#endif
#if defined (VECTORN_223) && !defined (SKIP_VECTORN_223)
    RunMe(223);
#endif
#if defined (VECTORN_224) && !defined (SKIP_VECTORN_224)
    RunMe(224);
#endif
#if defined (VECTORN_225) && !defined (SKIP_VECTORN_225)
    RunMe(225);
#endif
#if defined (VECTORN_226) && !defined (SKIP_VECTORN_226)
    RunMe(226);
#endif
#if defined (VECTORN_227) && !defined (SKIP_VECTORN_227)
    RunMe(227);
#endif
#if defined (VECTORN_228) && !defined (SKIP_VECTORN_228)
    RunMe(228);
#endif
#if defined (VECTORN_229) && !defined (SKIP_VECTORN_229)
    RunMe(229);
#endif
#if defined (VECTORN_230) && !defined (SKIP_VECTORN_230)
    RunMe(230);
#endif
#if defined (VECTORN_231) && !defined (SKIP_VECTORN_231)
    RunMe(231);
#endif
#if defined (VECTORN_232) && !defined (SKIP_VECTORN_232)
    RunMe(232);
#endif
#if defined (VECTORN_233) && !defined (SKIP_VECTORN_233)
    RunMe(233);
#endif
#if defined (VECTORN_234) && !defined (SKIP_VECTORN_234)
    RunMe(234);
#endif
#if defined (VECTORN_235) && !defined (SKIP_VECTORN_235)
    RunMe(235);
#endif
#if defined (VECTORN_236) && !defined (SKIP_VECTORN_236)
    RunMe(236);
#endif
#if defined (VECTORN_237) && !defined (SKIP_VECTORN_237)
    RunMe(237);
#endif
#if defined (VECTORN_238) && !defined (SKIP_VECTORN_238)
    RunMe(238);
#endif
#if defined (VECTORN_239) && !defined (SKIP_VECTORN_239)
    RunMe(239);
#endif
#if defined (VECTORN_240) && !defined (SKIP_VECTORN_240)
    RunMe(240);
#endif
#if defined (VECTORN_241) && !defined (SKIP_VECTORN_241)
    RunMe(241);
#endif
#if defined (VECTORN_242) && !defined (SKIP_VECTORN_242)
    RunMe(242);
#endif
#if defined (VECTORN_243) && !defined (SKIP_VECTORN_243)
    RunMe(243);
#endif
#if defined (VECTORN_244) && !defined (SKIP_VECTORN_244)
    RunMe(244);
#endif
#if defined (VECTORN_245) && !defined (SKIP_VECTORN_245)
    RunMe(245);
#endif
#if defined (VECTORN_246) && !defined (SKIP_VECTORN_246)
    RunMe(246);
#endif
#if defined (VECTORN_247) && !defined (SKIP_VECTORN_247)
    RunMe(247);
#endif
#if defined (VECTORN_248) && !defined (SKIP_VECTORN_248)
    RunMe(248);
#endif
#if defined (VECTORN_249) && !defined (SKIP_VECTORN_249)
    RunMe(249);
#endif
#if defined (VECTORN_250) && !defined (SKIP_VECTORN_250)
    RunMe(250);
#endif
#if defined (VECTORN_251) && !defined (SKIP_VECTORN_251)
    RunMe(251);
#endif
#if defined (VECTORN_252) && !defined (SKIP_VECTORN_252)
    RunMe(252);
#endif
#if defined (VECTORN_253) && !defined (SKIP_VECTORN_253)
    RunMe(253);
#endif
#if defined (VECTORN_254) && !defined (SKIP_VECTORN_254)
    RunMe(254);
#endif
#if defined (VECTORN_255) && !defined (SKIP_VECTORN_255)
    RunMe(255);
#endif

#undef VectorNMacro
#undef RunMe
#undef MacroArgs
