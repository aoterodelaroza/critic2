/*
Copyright (c) 2017 Alberto Otero de la Roza
<aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
<victor@fluor.quimica.uniovi.es>.

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* The spacegroup information comes from spacegroup.c, by Atsushi
Togo. See the license notice in src/spglib/spacegroup.c */
#include <stdlib.h>

#include "imgui.h"
#include "imgui_internal.h"
#include "guiapps.h"
#include "critic2.h"
#include "settings.h"
#include "draw.h"

#include "imguifilesystem.h"

// static function prototypes
static int spg_choose_menu(int mode);

// Variable definitions
bool structureinfo_window_h = false;
bool structurenew_window_h = false; 
int structureopen_window_h = 0; // 0 - hidden, 1 - molecule, 2 - crystal
bool console_window_h = false;

// space group types
static const struct {
  int number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_full[20];
  char international_short[11];
  char choice[6];
} spgtyp[] = {
{0,"","","","","",""},/*0*/
{1,"C1^1","P1","P1","P1","P1",""},/*1*/
{2,"Ci^1","-P1","P-1","P-1","P-1",""},/*2*/
{3,"C2^1","P2y","P2=P121","P121","P2","b"},/*3*/
{3,"C2^1","P2","P2=P112","P112","P2","c"},/*4*/
{3,"C2^1","P2x","P2=P211","P211","P2","a"},/*5*/
{4,"C2^2","P2yb","P2_1=P12_11","P12_11","P2_1","b"},/*6*/
{4,"C2^2","P2c","P2_1=P112_1","P112_1","P2_1","c"},/*7*/
{4,"C2^2","P2xa","P2_1=P2_111","P2_111","P2_1","a"},/*8*/
{5,"C2^3","C2y","C2=C121","C121","C2","b1"},/*9*/
{5,"C2^3","A2y","C2=A121","A121","C2","b2"},/*10*/
{5,"C2^3","I2y","C2=I121","I121","C2","b3"},/*11*/
{5,"C2^3","A2","C2=A112","A112","C2","c1"},/*12*/
{5,"C2^3","B2","C2=B112=B2","B112","C2","c2"},/*13*/
{5,"C2^3","I2","C2=I112","I112","C2","c3"},/*14*/
{5,"C2^3","B2x","C2=B211","B211","C2","a1"},/*15*/
{5,"C2^3","C2x","C2=C211","C211","C2","a2"},/*16*/
{5,"C2^3","I2x","C2=I211","I211","C2","a3"},/*17*/
{6,"Cs^1","P-2y","Pm=P1m1","P1m1","Pm","b"},/*18*/
{6,"Cs^1","P-2","Pm=P11m","P11m","Pm","c"},/*19*/
{6,"Cs^1","P-2x","Pm=Pm11","Pm11","Pm","a"},/*20*/
{7,"Cs^2","P-2yc","Pc=P1c1","P1c1","Pc","b1"},/*21*/
{7,"Cs^2","P-2yac","Pc=P1n1","P1n1","Pc","b2"},/*22*/
{7,"Cs^2","P-2ya","Pc=P1a1","P1a1","Pc","b3"},/*23*/
{7,"Cs^2","P-2a","Pc=P11a","P11a","Pc","c1"},/*24*/
{7,"Cs^2","P-2ab","Pc=P11n","P11n","Pc","c2"},/*25*/
{7,"Cs^2","P-2b","Pc=P11b=Pb","P11b","Pc","c3"},/*26*/
{7,"Cs^2","P-2xb","Pc=Pb11","Pb11","Pc","a1"},/*27*/
{7,"Cs^2","P-2xbc","Pc=Pn11","Pn11","Pc","a2"},/*28*/
{7,"Cs^2","P-2xc","Pc=Pc11","Pc11","Pc","a3"},/*29*/
{8,"Cs^3","C-2y","Cm=C1m1","C1m1","Cm","b1"},/*30*/
{8,"Cs^3","A-2y","Cm=A1m1","A1m1","Cm","b2"},/*31*/
{8,"Cs^3","I-2y","Cm=I1m1","I1m1","Cm","b3"},/*32*/
{8,"Cs^3","A-2","Cm=A11m","A11m","Cm","c1"},/*33*/
{8,"Cs^3","B-2","Cm=B11m=Bm","B11m","Cm","c2"},/*34*/
{8,"Cs^3","I-2","Cm=I11m","I11m","Cm","c3"},/*35*/
{8,"Cs^3","B-2x","Cm=Bm11","Bm11","Cm","a1"},/*36*/
{8,"Cs^3","C-2x","Cm=Cm11","Cm11","Cm","a2"},/*37*/
{8,"Cs^3","I-2x","Cm=Im11","Im11","Cm","a3"},/*38*/
{9,"Cs^4","C-2yc","Cc=C1c1","C1c1","Cc","b1"},/*39*/
{9,"Cs^4","A-2yac","Cc=A1n1","A1n1","Cc","b2"},/*40*/
{9,"Cs^4","I-2ya","Cc=I1a1","I1a1","Cc","b3"},/*41*/
{9,"Cs^4","A-2ya","Cc=A1a1","A1a1","Cc","-b1"},/*42*/
{9,"Cs^4","C-2ybc","Cc=C1n1","C1n1","Cc","-b2"},/*43*/
{9,"Cs^4","I-2yc","Cc=I1c1","I1c1","Cc","-b3"},/*44*/
{9,"Cs^4","A-2a","Cc=A11a","A11a","Cc","c1"},/*45*/
{9,"Cs^4","B-2bc","Cc=B11n","B11n","Cc","c2"},/*46*/
{9,"Cs^4","I-2b","Cc=I11b","I11b","Cc","c3"},/*47*/
{9,"Cs^4","B-2b","Cc=B11b=Bb","B11b","Cc","-c1"},/*48*/
{9,"Cs^4","A-2ac","Cc=A11n","A11n","Cc","-c2"},/*49*/
{9,"Cs^4","I-2a","Cc=I11a","I11a","Cc","-c3"},/*50*/
{9,"Cs^4","B-2xb","Cc=Bb11","Bb11","Cc","a1"},/*51*/
{9,"Cs^4","C-2xbc","Cc=Cn11","Cn11","Cc","a2"},/*52*/
{9,"Cs^4","I-2xc","Cc=Ic11","Ic11","Cc","a3"},/*53*/
{9,"Cs^4","C-2xc","Cc=Cc11","Cc11","Cc","-a1"},/*54*/
{9,"Cs^4","B-2xbc","Cc=Bn11","Bn11","Cc","-a2"},/*55*/
{9,"Cs^4","I-2xb","Cc=Ib11","Ib11","Cc","-a3"},/*56*/
{10,"C2h^1","-P2y","P2/m=P12/m1","P12/m1","P2/m","b"},/*57*/
{10,"C2h^1","-P2","P2/m=P112/m","P112/m","P2/m","c"},/*58*/
{10,"C2h^1","-P2x","P2/m=P2/m11","P2/m11","P2/m","a"},/*59*/
{11,"C2h^2","-P2yb","P2_1/m=P12_1/m1","P12_1/m1","P2_1/m","b"},/*60*/
{11,"C2h^2","-P2c","P2_1/m=P112_1/m","P112_1/m","P2_1/m","c"},/*61*/
{11,"C2h^2","-P2xa","P2_1/m=P2_1/m11","P2_1/m11","P2_1/m","a"},/*62*/
{12,"C2h^3","-C2y","C2/m=C12/m1","C12/m1","C2/m","b1"},/*63*/
{12,"C2h^3","-A2y","C2/m=A12/m1","A12/m1","C2/m","b2"},/*64*/
{12,"C2h^3","-I2y","C2/m=I12/m1","I12/m1","C2/m","b3"},/*65*/
{12,"C2h^3","-A2","C2/m=A112/m","A112/m","C2/m","c1"},/*66*/
{12,"C2h^3","-B2","C2/m=B112/m=B2/m","B112/m","C2/m","c2"},/*67*/
{12,"C2h^3","-I2","C2/m=I112/m","I112/m","C2/m","c3"},/*68*/
{12,"C2h^3","-B2x","C2/m=B2/m11","B2/m11","C2/m","a1"},/*69*/
{12,"C2h^3","-C2x","C2/m=C2/m11","C2/m11","C2/m","a2"},/*70*/
{12,"C2h^3","-I2x","C2/m=I2/m11","I2/m11","C2/m","a3"},/*71*/
{13,"C2h^4","-P2yc","P2/c=P12/c1","P12/c1","P2/c","b1"},/*72*/
{13,"C2h^4","-P2yac","P2/c=P12/n1","P12/n1","P2/c","b2"},/*73*/
{13,"C2h^4","-P2ya","P2/c=P12/a1","P12/a1","P2/c","b3"},/*74*/
{13,"C2h^4","-P2a","P2/c=P112/a","P112/a","P2/c","c1"},/*75*/
{13,"C2h^4","-P2ab","P2/c=P112/n","P112/n","P2/c","c2"},/*76*/
{13,"C2h^4","-P2b","P2/c=P112/b=P2/b","P112/b","P2/c","c3"},/*77*/
{13,"C2h^4","-P2xb","P2/c=P2/b11","P2/b11","P2/c","a1"},/*78*/
{13,"C2h^4","-P2xbc","P2/c=P2/n11","P2/n11","P2/c","a2"},/*79*/
{13,"C2h^4","-P2xc","P2/c=P2/c11","P2/c11","P2/c","a3"},/*80*/
{14,"C2h^5","-P2ybc","P2_1/c=P12_1/c1","P12_1/c1","P2_1/c","b1"},/*81*/
{14,"C2h^5","-P2yn","P2_1/c=P12_1/n1","P12_1/n1","P2_1/c","b2"},/*82*/
{14,"C2h^5","-P2yab","P2_1/c=P12_1/a1","P12_1/a1","P2_1/c","b3"},/*83*/
{14,"C2h^5","-P2ac","P2_1/c=P112_1/a","P112_1/a","P2_1/c","c1"},/*84*/
{14,"C2h^5","-P2n","P2_1/c=P112_1/n","P112_1/n","P2_1/c","c2"},/*85*/
{14,"C2h^5","-P2bc","P2_1/c=P112_1/b=P2_1/b","P112_1/b","P2_1/c","c3"},/*86*/
{14,"C2h^5","-P2xab","P2_1/c=P2_1/b11","P2_1/b11","P2_1/c","a1"},/*87*/
{14,"C2h^5","-P2xn","P2_1/c=P2_1/n11","P2_1/n11","P2_1/c","a2"},/*88*/
{14,"C2h^5","-P2xac","P2_1/c=P2_1/c11","P2_1/c11","P2_1/c","a3"},/*89*/
{15,"C2h^6","-C2yc","C2/c=C12/c1","C12/c1","C2/c","b1"},/*90*/
{15,"C2h^6","-A2yac","C2/c=A12/n1","A12/n1","C2/c","b2"},/*91*/
{15,"C2h^6","-I2ya","C2/c=I12/a1","I12/a1","C2/c","b3"},/*92*/
{15,"C2h^6","-A2ya","C2/c=A12/a1","A12/a1","C2/c","-b1"},/*93*/
{15,"C2h^6","-C2ybc","C2/c=C12/n1","C12/n1","C2/c","-b2"},/*94*/
{15,"C2h^6","-I2yc","C2/c=I12/c1","I12/c1","C2/c","-b3"},/*95*/
{15,"C2h^6","-A2a","C2/c=A112/a","A112/a","C2/c","c1"},/*96*/
{15,"C2h^6","-B2bc","C2/c=B112/n","B112/n","C2/c","c2"},/*97*/
{15,"C2h^6","-I2b","C2/c=I112/b","I112/b","C2/c","c3"},/*98*/
{15,"C2h^6","-B2b","C2/c=B112/b=B2/b","B112/b","C2/c","-c1"},/*99*/
{15,"C2h^6","-A2ac","C2/c=A112/n","A112/n","C2/c","-c2"},/*100*/
{15,"C2h^6","-I2a","C2/c=I112/a","I112/a","C2/c","-c3"},/*101*/
{15,"C2h^6","-B2xb","C2/c=B2/b11","B2/b11","C2/c","a1"},/*102*/
{15,"C2h^6","-C2xbc","C2/c=C2/n11","C2/n11","C2/c","a2"},/*103*/
{15,"C2h^6","-I2xc","C2/c=I2/c11","I2/c11","C2/c","a3"},/*104*/
{15,"C2h^6","-C2xc","C2/c=C2/c11","C2/c11","C2/c","-a1"},/*105*/
{15,"C2h^6","-B2xbc","C2/c=B2/n11","B2/n11","C2/c","-a2"},/*106*/
{15,"C2h^6","-I2xb","C2/c=I2/b11","I2/b11","C2/c","-a3"},/*107*/
{16,"D2^1","P22","P222","P222","P222",""},/*108*/
{17,"D2^2","P2c2","P222_1","P222_1","P222_1",""},/*109*/
{17,"D2^2","P2a2a","P2_122","P2_122","P2_122","cab"},/*110*/
{17,"D2^2","P22b","P22_12","P22_12","P22_12","bca"},/*111*/
{18,"D2^3","P22ab","P2_12_12","P2_12_12","P2_12_12",""},/*112*/
{18,"D2^3","P2bc2","P22_12_1","P22_12_1","P22_12_1","cab"},/*113*/
{18,"D2^3","P2ac2ac","P2_122_1","P2_122_1","P2_122_1","bca"},/*114*/
{19,"D2^4","P2ac2ab","P2_12_12_1","P2_12_12_1","P2_12_12_1",""},/*115*/
{20,"D2^5","C2c2","C222_1","C222_1","C222_1",""},/*116*/
{20,"D2^5","A2a2a","A2_122","A2_122","A2_122","cab"},/*117*/
{20,"D2^5","B22b","B22_12","B22_12","B22_12","bca"},/*118*/
{21,"D2^6","C22","C222","C222","C222",""},/*119*/
{21,"D2^6","A22","A222","A222","A222","cab"},/*120*/
{21,"D2^6","B22","B222","B222","B222","bca"},/*121*/
{22,"D2^7","F22","F222","F222","F222",""},/*122*/
{23,"D2^8","I22","I222","I222","I222",""},/*123*/
{24,"D2^9","I2b2c","I2_12_12_1","I2_12_12_1","I2_12_12_1",""},/*124*/
{25,"C2v^1","P2-2","Pmm2","Pmm2","Pmm2",""},/*125*/
{25,"C2v^1","P-22","P2mm","P2mm","P2mm","cab"},/*126*/
{25,"C2v^1","P-2-2","Pm2m","Pm2m","Pm2m","bca"},/*127*/
{26,"C2v^2","P2c-2","Pmc2_1","Pmc2_1","Pmc2_1",""},/*128*/
{26,"C2v^2","P2c-2c","Pcm2_1","Pcm2_1","Pcm2_1","ba-c"},/*129*/
{26,"C2v^2","P-2a2a","P2_1ma","P2_1ma","P2_1ma","cab"},/*130*/
{26,"C2v^2","P-22a","P2_1am","P2_1am","P2_1am","-cba"},/*131*/
{26,"C2v^2","P-2-2b","Pb2_1m","Pb2_1m","Pb2_1m","bca"},/*132*/
{26,"C2v^2","P-2b-2","Pm2_1b","Pm2_1b","Pm2_1b","a-cb"},/*133*/
{27,"C2v^3","P2-2c","Pcc2","Pcc2","Pcc2",""},/*134*/
{27,"C2v^3","P-2a2","P2aa","P2aa","P2aa","cab"},/*135*/
{27,"C2v^3","P-2b-2b","Pb2b","Pb2b","Pb2b","bca"},/*136*/
{28,"C2v^4","P2-2a","Pma2","Pma2","Pma2",""},/*137*/
{28,"C2v^4","P2-2b","Pbm2","Pbm2","Pbm2","ba-c"},/*138*/
{28,"C2v^4","P-2b2","P2mb","P2mb","P2mb","cab"},/*139*/
{28,"C2v^4","P-2c2","P2cm","P2cm","P2cm","-cba"},/*140*/
{28,"C2v^4","P-2c-2c","Pc2m","Pc2m","Pc2m","bca"},/*141*/
{28,"C2v^4","P-2a-2a","Pm2a","Pm2a","Pm2a","a-cb"},/*142*/
{29,"C2v^5","P2c-2ac","Pca2_1","Pca2_1","Pca2_1",""},/*143*/
{29,"C2v^5","P2c-2b","Pbc2_1","Pbc2_1","Pbc2_1","ba-c"},/*144*/
{29,"C2v^5","P-2b2a","P2_1ab","P2_1ab","P2_1ab","cab"},/*145*/
{29,"C2v^5","P-2ac2a","P2_1ca","P2_1ca","P2_1ca","-cba"},/*146*/
{29,"C2v^5","P-2bc-2c","Pc2_1b","Pc2_1b","Pc2_1b","bca"},/*147*/
{29,"C2v^5","P-2a-2ab","Pb2_1a","Pb2_1a","Pb2_1a","a-cb"},/*148*/
{30,"C2v^6","P2-2bc","Pnc2","Pnc2","Pnc2",""},/*149*/
{30,"C2v^6","P2-2ac","Pcn2","Pcn2","Pcn2","ba-c"},/*150*/
{30,"C2v^6","P-2ac2","P2na","P2na","P2na","cab"},/*151*/
{30,"C2v^6","P-2ab2","P2an","P2an","P2an","-cba"},/*152*/
{30,"C2v^6","P-2ab-2ab","Pb2n","Pb2n","Pb2n","bca"},/*153*/
{30,"C2v^6","P-2bc-2bc","Pn2b","Pn2b","Pn2b","a-cb"},/*154*/
{31,"C2v^7","P2ac-2","Pmn2_1","Pmn2_1","Pmn2_1",""},/*155*/
{31,"C2v^7","P2bc-2bc","Pnm2_1","Pnm2_1","Pnm2_1","ba-c"},/*156*/
{31,"C2v^7","P-2ab2ab","P2_1mn","P2_1mn","P2_1mn","cab"},/*157*/
{31,"C2v^7","P-22ac","P2_1nm","P2_1nm","P2_1nm","-cba"},/*158*/
{31,"C2v^7","P-2-2bc","Pn2_1m","Pn2_1m","Pn2_1m","bca"},/*159*/
{31,"C2v^7","P-2ab-2","Pm2_1n","Pm2_1n","Pm2_1n","a-cb"},/*160*/
{32,"C2v^8","P2-2ab","Pba2","Pba2","Pba2",""},/*161*/
{32,"C2v^8","P-2bc2","P2cb","P2cb","P2cb","cab"},/*162*/
{32,"C2v^8","P-2ac-2ac","Pc2a","Pc2a","Pc2a","bca"},/*163*/
{33,"C2v^9","P2c-2n","Pna2_1","Pna2_1","Pna2_1",""},/*164*/
{33,"C2v^9","P2c-2ab","Pbn2_1","Pbn2_1","Pbn2_1","ba-c"},/*165*/
{33,"C2v^9","P-2bc2a","P2_1nb","P2_1nb","P2_1nb","cab"},/*166*/
{33,"C2v^9","P-2n2a","P2_1cn","P2_1cn","P2_1cn","-cba"},/*167*/
{33,"C2v^9","P-2n-2ac","Pc2_1n","Pc2_1n","Pc2_1n","bca"},/*168*/
{33,"C2v^9","P-2ac-2n","Pn2_1a","Pn2_1a","Pn2_1a","a-cb"},/*169*/
{34,"C2v^10","P2-2n","Pnn2","Pnn2","Pnn2",""},/*170*/
{34,"C2v^10","P-2n2","P2nn","P2nn","P2nn","cab"},/*171*/
{34,"C2v^10","P-2n-2n","Pn2n","Pn2n","Pn2n","bca"},/*172*/
{35,"C2v^11","C2-2","Cmm2","Cmm2","Cmm2",""},/*173*/
{35,"C2v^11","A-22","A2mm","A2mm","A2mm","cab"},/*174*/
{35,"C2v^11","B-2-2","Bm2m","Bm2m","Bm2m","bca"},/*175*/
{36,"C2v^12","C2c-2","Cmc2_1","Cmc2_1","Cmc2_1",""},/*176*/
{36,"C2v^12","C2c-2c","Ccm2_1","Ccm2_1","Ccm2_1","ba-c"},/*177*/
{36,"C2v^12","A-2a2a","A2_1ma","A2_1ma","A2_1ma","cab"},/*178*/
{36,"C2v^12","A-22a","A2_1am","A2_1am","A2_1am","-cba"},/*179*/
{36,"C2v^12","B-2-2b","Bb2_1m","Bb2_1m","Bb2_1m","bca"},/*180*/
{36,"C2v^12","B-2b-2","Bm2_1b","Bm2_1b","Bm2_1b","a-cb"},/*181*/
{37,"C2v^13","C2-2c","Ccc2","Ccc2","Ccc2",""},/*182*/
{37,"C2v^13","A-2a2","A2aa","A2aa","A2aa","cab"},/*183*/
{37,"C2v^13","B-2b-2b","Bb2b","Bb2b","Bb2b","bca"},/*184*/
{38,"C2v^14","A2-2","Amm2","Amm2","Amm2",""},/*185*/
{38,"C2v^14","B2-2","Bmm2","Bmm2","Bmm2","ba-c"},/*186*/
{38,"C2v^14","B-22","B2mm","B2mm","B2mm","cab"},/*187*/
{38,"C2v^14","C-22","C2mm","C2mm","C2mm","-cba"},/*188*/
{38,"C2v^14","C-2-2","Cm2m","Cm2m","Cm2m","bca"},/*189*/
{38,"C2v^14","A-2-2","Am2m","Am2m","Am2m","a-cb"},/*190*/
{39,"C2v^15","A2-2c","Aem2","Aem2","Aem2",""},/*191*/
{39,"C2v^15","B2-2c","Bme2","Bme2","Bme2","ba-c"},/*192*/
{39,"C2v^15","B-2c2","B2em","B2em","B2em","cab"},/*193*/
{39,"C2v^15","C-2b2","C2me","C2me","C2me","-cba"},/*194*/
{39,"C2v^15","C-2b-2b","Cm2e","Cm2e","Cm2e","bca"},/*195*/
{39,"C2v^15","A-2c-2c","Ae2m","Ae2m","Ae2m","a-cb"},/*196*/
{40,"C2v^16","A2-2a","Ama2","Ama2","Ama2",""},/*197*/
{40,"C2v^16","B2-2b","Bbm2","Bbm2","Bbm2","ba-c"},/*198*/
{40,"C2v^16","B-2b2","B2mb","B2mb","B2mb","cab"},/*199*/
{40,"C2v^16","C-2c2","C2cm","C2cm","C2cm","-cba"},/*200*/
{40,"C2v^16","C-2c-2c","Cc2m","Cc2m","Cc2m","bca"},/*201*/
{40,"C2v^16","A-2a-2a","Am2a","Am2a","Am2a","a-cb"},/*202*/
{41,"C2v^17","A2-2ac","Aea2","Aea2","Aea2",""},/*203*/
{41,"C2v^17","B2-2bc","Bbe2","Bbe2","Bbe2","ba-c"},/*204*/
{41,"C2v^17","B-2bc2","B2eb","B2eb","B2eb","cab"},/*205*/
{41,"C2v^17","C-2bc2","C2ce","C2ce","C2ce","-cba"},/*206*/
{41,"C2v^17","C-2bc-2bc","Cc2e","Cc2e","Cc2e","bca"},/*207*/
{41,"C2v^17","A-2ac-2ac","Ae2a","Ae2a","Ae2a","a-cb"},/*208*/
{42,"C2v^18","F2-2","Fmm2","Fmm2","Fmm2",""},/*209*/
{42,"C2v^18","F-22","F2mm","F2mm","F2mm","cab"},/*210*/
{42,"C2v^18","F-2-2","Fm2m","Fm2m","Fm2m","bca"},/*211*/
{43,"C2v^19","F2-2d","Fdd2","Fdd2","Fdd2",""},/*212*/
{43,"C2v^19","F-2d2","F2dd","F2dd","F2dd","cab"},/*213*/
{43,"C2v^19","F-2d-2d","Fd2d","Fd2d","Fd2d","bca"},/*214*/
{44,"C2v^20","I2-2","Imm2","Imm2","Imm2",""},/*215*/
{44,"C2v^20","I-22","I2mm","I2mm","I2mm","cab"},/*216*/
{44,"C2v^20","I-2-2","Im2m","Im2m","Im2m","bca"},/*217*/
{45,"C2v^21","I2-2c","Iba2","Iba2","Iba2",""},/*218*/
{45,"C2v^21","I-2a2","I2cb","I2cb","I2cb","cab"},/*219*/
{45,"C2v^21","I-2b-2b","Ic2a","Ic2a","Ic2a","bca"},/*220*/
{46,"C2v^22","I2-2a","Ima2","Ima2","Ima2",""},/*221*/
{46,"C2v^22","I2-2b","Ibm2","Ibm2","Ibm2","ba-c"},/*222*/
{46,"C2v^22","I-2b2","I2mb","I2mb","I2mb","cab"},/*223*/
{46,"C2v^22","I-2c2","I2cm","I2cm","I2cm","-cba"},/*224*/
{46,"C2v^22","I-2c-2c","Ic2m","Ic2m","Ic2m","bca"},/*225*/
{46,"C2v^22","I-2a-2a","Im2a","Im2a","Im2a","a-cb"},/*226*/
{47,"D2h^1","-P22","Pmmm","P2/m2/m2/m","Pmmm",""},/*227*/
{48,"D2h^2","P22-1n","Pnnn","P2/n2/n2/n","Pnnn","1"},/*228*/
{48,"D2h^2","-P2ab2bc","Pnnn","P2/n2/n2/n","Pnnn","2"},/*229*/
{49,"D2h^3","-P22c","Pccm","P2/c2/c2/m","Pccm",""},/*230*/
{49,"D2h^3","-P2a2","Pmaa","P2/m2/a2/a","Pmaa","cab"},/*231*/
{49,"D2h^3","-P2b2b","Pbmb","P2/b2/m2/b","Pbmb","bca"},/*232*/
{50,"D2h^4","P22-1ab","Pban","P2/b2/a2/n","Pban","1"},/*233*/
{50,"D2h^4","-P2ab2b","Pban","P2/b2/a2/n","Pban","2"},/*234*/
{50,"D2h^4","P22-1bc","Pncb","P2/n2/c2/b","Pncb","1cab"},/*235*/
{50,"D2h^4","-P2b2bc","Pncb","P2/n2/c2/b","Pncb","2cab"},/*236*/
{50,"D2h^4","P22-1ac","Pcna","P2/c2/n2/a","Pcna","1bca"},/*237*/
{50,"D2h^4","-P2a2c","Pcna","P2/c2/n2/a","Pcna","2bca"},/*238*/
{51,"D2h^5","-P2a2a","Pmma","P2_1/m2/m2/a","Pmma",""},/*239*/
{51,"D2h^5","-P2b2","Pmmb","P2/m2_1/m2/b","Pmmb","ba-c"},/*240*/
{51,"D2h^5","-P22b","Pbmm","P2/b2_1/m2/m","Pbmm","cab"},/*241*/
{51,"D2h^5","-P2c2c","Pcmm","P2/c2/m2_1/m","Pcmm","-cba"},/*242*/
{51,"D2h^5","-P2c2","Pmcm","P2/m2/c2_1/m","Pmcm","bca"},/*243*/
{51,"D2h^5","-P22a","Pmam","P2_1/m2/a2/m","Pmam","a-cb"},/*244*/
{52,"D2h^6","-P2a2bc","Pnna","P2/n2_1/n2/a","Pnna",""},/*245*/
{52,"D2h^6","-P2b2n","Pnnb","P2_1/n2/n2/b","Pnnb","ba-c"},/*246*/
{52,"D2h^6","-P2n2b","Pbnn","P2/b2/n2_1/n","Pbnn","cab"},/*247*/
{52,"D2h^6","-P2ab2c","Pcnn","P2/c2_1/n2/n","Pcnn","-cba"},/*248*/
{52,"D2h^6","-P2ab2n","Pncn","P2_1/n2/c2/n","Pncn","bca"},/*249*/
{52,"D2h^6","-P2n2bc","Pnan","P2/n2/a2_1/n","Pnan","a-cb"},/*250*/
{53,"D2h^7","-P2ac2","Pmna","P2/m2/n2_1/a","Pmna",""},/*251*/
{53,"D2h^7","-P2bc2bc","Pnmb","P2/n2/m2_1/b","Pnmb","ba-c"},/*252*/
{53,"D2h^7","-P2ab2ab","Pbmn","P2_1/b2/m2/n","Pbmn","cab"},/*253*/
{53,"D2h^7","-P22ac","Pcnm","P2_1/c2/n2/m","Pcnm","-cba"},/*254*/
{53,"D2h^7","-P22bc","Pncm","P2/n2_1/c2/m","Pncm","bca"},/*255*/
{53,"D2h^7","-P2ab2","Pman","P2/m2_1/a2/n","Pman","a-cb"},/*256*/
{54,"D2h^8","-P2a2ac","Pcca","P2_1/c2/c2/a","Pcca",""},/*257*/
{54,"D2h^8","-P2b2c","Pccb","P2/c2_1/c2/b","Pccb","ba-c"},/*258*/
{54,"D2h^8","-P2a2b","Pbaa","P2/b2_1/a2/a","Pbaa","cab"},/*259*/
{54,"D2h^8","-P2ac2c","Pcaa","P2/c2/a2_1/a","Pcaa","-cba"},/*260*/
{54,"D2h^8","-P2bc2b","Pbcb","P2/b2/c2_1/b","Pbcb","bca"},/*261*/
{54,"D2h^8","-P2b2ab","Pbab","P2_1/b2/a2/b","Pbab","a-cb"},/*262*/
{55,"D2h^9","-P22ab","Pbam","P2_1/b2_1/a2/m","Pbam",""},/*263*/
{55,"D2h^9","-P2bc2","Pmcb","P2/m2_1/c2_1/b","Pmcb","cab"},/*264*/
{55,"D2h^9","-P2ac2ac","Pcma","P2_1/c2/m2_1/a","Pcma","bca"},/*265*/
{56,"D2h^10","-P2ab2ac","Pccn","P2_1/c2_1/c2/n","Pccn",""},/*266*/
{56,"D2h^10","-P2ac2bc","Pnaa","P2/n2_1/a2_1/a","Pnaa","cab"},/*267*/
{56,"D2h^10","-P2bc2ab","Pbnb","P2_1/b2/n2_1/b","Pbnb","bca"},/*268*/
{57,"D2h^11","-P2c2b","Pbcm","P2/b2_1/c2_1/m","Pbcm",""},/*269*/
{57,"D2h^11","-P2c2ac","Pcam","P2_1/c2/a2_1/m","Pcam","ba-c"},/*270*/
{57,"D2h^11","-P2ac2a","Pmca","P2_1/m2/c2_1/a","Pmca","cab"},/*271*/
{57,"D2h^11","-P2b2a","Pmab","P2_1/m2_1/a2/b","Pmab","-cba"},/*272*/
{57,"D2h^11","-P2a2ab","Pbma","P2_1/b2_1/m2/a","Pbma","bca"},/*273*/
{57,"D2h^11","-P2bc2c","Pcmb","P2/c2_1/m2_1/b","Pcmb","a-cb"},/*274*/
{58,"D2h^12","-P22n","Pnnm","P2_1/n2_1/n2/m","Pnnm",""},/*275*/
{58,"D2h^12","-P2n2","Pmnn","P2/m2_1/n2_1/n","Pmnn","cab"},/*276*/
{58,"D2h^12","-P2n2n","Pnmn","P2_1/n2/m2_1/n","Pnmn","bca"},/*277*/
{59,"D2h^13","P22ab-1ab","Pmmn","P2_1/m2_1/m2/n","Pmmn","1"},/*278*/
{59,"D2h^13","-P2ab2a","Pmmn","P2_1/m2_1/m2/n","Pmmn","2"},/*279*/
{59,"D2h^13","P2bc2-1bc","Pnmm","P2/n2_1/m2_1/m","Pnmm","1cab"},/*280*/
{59,"D2h^13","-P2c2bc","Pnmm","P2/n2_1/m2_1/m","Pnmm","2cab"},/*281*/
{59,"D2h^13","P2ac2ac-1ac","Pmnm","P2_1/m2/n2_1/m","Pmnm","1bca"},/*282*/
{59,"D2h^13","-P2c2a","Pmnm","P2_1/m2/n2_1/m","Pmnm","2bca"},/*283*/
{60,"D2h^14","-P2n2ab","Pbcn","P2_1/b2/c2_1/n","Pbcn",""},/*284*/
{60,"D2h^14","-P2n2c","Pcan","P2/c2_1/a2_1/n","Pcan","ba-c"},/*285*/
{60,"D2h^14","-P2a2n","Pnca","P2_1/n2_1/c2/a","Pnca","cab"},/*286*/
{60,"D2h^14","-P2bc2n","Pnab","P2_1/n2/a2_1/b","Pnab","-cba"},/*287*/
{60,"D2h^14","-P2ac2b","Pbna","P2/b2_1/n2_1/a","Pbna","bca"},/*288*/
{60,"D2h^14","-P2b2ac","Pcnb","P2_1/c2_1/n2/b","Pcnb","a-cb"},/*289*/
{61,"D2h^15","-P2ac2ab","Pbca","P2_1/b2_1/c2_1/a","Pbca",""},/*290*/
{61,"D2h^15","-P2bc2ac","Pcab","P2_1/c2_1/a2_1/b","Pcab","ba-c"},/*291*/
{62,"D2h^16","-P2ac2n","Pnma","P2_1/n2_1/m2_1/a","Pnma",""},/*292*/
{62,"D2h^16","-P2bc2a","Pmnb","P2_1/m2_1/n2_1/b","Pmnb","ba-c"},/*293*/
{62,"D2h^16","-P2c2ab","Pbnm","P2_1/b2_1/n2_1/m","Pbnm","cab"},/*294*/
{62,"D2h^16","-P2n2ac","Pcmn","P2_1/c2_1/m2_1/n","Pcmn","-cba"},/*295*/
{62,"D2h^16","-P2n2a","Pmcn","P2_1/m2_1/c2_1/n","Pmcn","bca"},/*296*/
{62,"D2h^16","-P2c2n","Pnam","P2_1/n2_1/a2_1/m","Pnam","a-cb"},/*297*/
{63,"D2h^17","-C2c2","Cmcm","C2/m2/c2_1/m","Cmcm",""},/*298*/
{63,"D2h^17","-C2c2c","Ccmm","C2/c2/m2_1/m","Ccmm","ba-c"},/*299*/
{63,"D2h^17","-A2a2a","Amma","A2_1/m2/m2/a","Amma","cab"},/*300*/
{63,"D2h^17","-A22a","Amam","A2_1/m2/a2/m","Amam","-cba"},/*301*/
{63,"D2h^17","-B22b","Bbmm","B2/b2_1/m2/m","Bbmm","bca"},/*302*/
{63,"D2h^17","-B2b2","Bmmb","B2/m2_1/m2/b","Bmmb","a-cb"},/*303*/
{64,"D2h^18","-C2bc2","Cmce","C2/m2/c2_1/e","Cmce",""},/*304*/
{64,"D2h^18","-C2bc2bc","Ccme","C2/c2/m2_1/e","Ccme","ba-c"},/*305*/
{64,"D2h^18","-A2ac2ac","Aema","A2_1/e2/m2/a","Aema","cab"},/*306*/
{64,"D2h^18","-A22ac","Aeam","A2_1/e2/a2/m","Aeam","-cba"},/*307*/
{64,"D2h^18","-B22bc","Bbem","B2/b2_1/e2/m","Bbem","bca"},/*308*/
{64,"D2h^18","-B2bc2","Bmeb","B2/m2_1/e2/b","Bmeb","a-cb"},/*309*/
{65,"D2h^19","-C22","Cmmm","C2/m2/m2/m","Cmmm",""},/*310*/
{65,"D2h^19","-A22","Ammm","A2/m2/m2/m","Ammm","cab"},/*311*/
{65,"D2h^19","-B22","Bmmm","B2/m2/m2/m","Bmmm","bca"},/*312*/
{66,"D2h^20","-C22c","Cccm","C2/c2/c2/m","Cccm",""},/*313*/
{66,"D2h^20","-A2a2","Amaa","A2/m2/a2/a","Amaa","cab"},/*314*/
{66,"D2h^20","-B2b2b","Bbmb","B2/b2/m2/b","Bbmb","bca"},/*315*/
{67,"D2h^21","-C2b2","Cmme","C2/m2/m2/e","Cmme",""},/*316*/
{67,"D2h^21","-C2b2b","Cmme","C2/m2/m2/e","Cmme","ba-c"},/*317*/
{67,"D2h^21","-A2c2c","Aemm","A2/e2/m2/m","Aemm","cab"},/*318*/
{67,"D2h^21","-A22c","Aemm","A2/e2/m2/m","Aemm","-cba"},/*319*/
{67,"D2h^21","-B22c","Bmem","B2/m2/e2/m","Bmem","bca"},/*320*/
{67,"D2h^21","-B2c2","Bmem","B2/m2/e2/m","Bmem","a-cb"},/*321*/
{68,"D2h^22","C22-1bc","Ccce","C2/c2/c2/e","Ccce","1"},/*322*/
{68,"D2h^22","-C2b2bc","Ccce","C2/c2/c2/e","Ccce","2"},/*323*/
{68,"D2h^22","C22-1bc","Ccce","C2/c2/c2/e","Ccce","1ba"},/*324*/
{68,"D2h^22","-C2b2c","Ccce","C2/c2/c2/e","Ccce","2ba"},/*325*/
{68,"D2h^22","A22-1ac","Aeaa","A2/e2/a2/a","Aeaa","1cab"},/*326*/
{68,"D2h^22","-A2a2c","Aeaa","A2/e2/a2/a","Aeaa","2cab"},/*327*/
{68,"D2h^22","A22-1ac","Aeaa","A2/e2/a2/a","Aeaa","1-cba"},/*328*/
{68,"D2h^22","-A2ac2c","Aeaa","A2/e2/a2/a","Aeaa","2-cba"},/*329*/
{68,"D2h^22","B22-1bc","Bbeb","B2/b2/e2/b","Bbeb","1bca"},/*330*/
{68,"D2h^22","-B2bc2b","Bbcb","B2/b2/e2/b","Bbcb","2bca"},/*331*/
{68,"D2h^22","B22-1bc","Bbeb","B2/b2/e2/b","Bbeb","1a-cb"},/*332*/
{68,"D2h^22","-B2b2bc","Bbeb","B2/b2/e2/b","Bbeb","2a-cb"},/*333*/
{69,"D2h^23","-F22","Fmmm","F2/m2/m2/m","Fmmm",""},/*334*/
{70,"D2h^24","F22-1d","Fddd","F2/d2/d2/d","Fddd","1"},/*335*/
{70,"D2h^24","-F2uv2vw","Fddd","F2/d2/d2/d","Fddd","2"},/*336*/
{71,"D2h^25","-I22","Immm","I2/m2/m2/m","Immm",""},/*337*/
{72,"D2h^26","-I22c","Ibam","I2/b2/a2/m","Ibam",""},/*338*/
{72,"D2h^26","-I2a2","Imcb","I2/m2/c2/b","Imcb","cab"},/*339*/
{72,"D2h^26","-I2b2b","Icma","I2/c2/m2/a","Icma","bca"},/*340*/
{73,"D2h^27","-I2b2c","Ibca","I2/b2/c2/a","Ibca",""},/*341*/
{73,"D2h^27","-I2a2b","Icab","I2/c2/a2/b","Icab","ba-c"},/*342*/
{74,"D2h^28","-I2b2","Imma","I2/m2/m2/a","Imma",""},/*343*/
{74,"D2h^28","-I2a2a","Immb","I2/m2/m2/b","Immb","ba-c"},/*344*/
{74,"D2h^28","-I2c2c","Ibmm","I2/b2/m2/m","Ibmm","cab"},/*345*/
{74,"D2h^28","-I22b","Icmm","I2/c2/m2/m","Icmm","-cba"},/*346*/
{74,"D2h^28","-I22a","Imcm","I2/m2/c2/m","Imcm","bca"},/*347*/
{74,"D2h^28","-I2c2","Imam","I2/m2/a2/m","Imam","a-cb"},/*348*/
{75,"C4^1","P4","P4","P4","P4",""},/*349*/
{76,"C4^2","P4w","P4_1","P4_1","P4_1",""},/*350*/
{77,"C4^3","P4c","P4_2","P4_2","P4_2",""},/*351*/
{78,"C4^4","P4cw","P4_3","P4_3","P4_3",""},/*352*/
{79,"C4^5","I4","I4","I4","I4",""},/*353*/
{80,"C4^6","I4bw","I4_1","I4_1","I4_1",""},/*354*/
{81,"S4^1","P-4","P-4","P-4","P-4",""},/*355*/
{82,"S4^2","I-4","I-4","I-4","I-4",""},/*356*/
{83,"C4h^1","-P4","P4/m","P4/m","P4/m",""},/*357*/
{84,"C4h^2","-P4c","P4_2/m","P4_2/m","P4_2/m",""},/*358*/
{85,"C4h^3","P4ab-1ab","P4/n","P4/n","P4/n","1"},/*359*/
{85,"C4h^3","-P4a","P4/n","P4/n","P4/n","2"},/*360*/
{86,"C4h^4","P4n-1n","P4_2/n","P4_2/n","P4_2/n","1"},/*361*/
{86,"C4h^4","-P4bc","P4_2/n","P4_2/n","P4_2/n","2"},/*362*/
{87,"C4h^5","-I4","I4/m","I4/m","I4/m",""},/*363*/
{88,"C4h^6","I4bw-1bw","I4_1/a","I4_1/a","I4_1/a","1"},/*364*/
{88,"C4h^6","-I4ad","I4_1/a","I4_1/a","I4_1/a","2"},/*365*/
{89,"D4^1","P42","P422","P422","P422",""},/*366*/
{90,"D4^2","P4ab2ab","P42_12","P42_12","P42_12",""},/*367*/
{91,"D4^3","P4w2c","P4_122","P4_122","P4_122",""},/*368*/
{92,"D4^4","P4abw2nw","P4_12_12","P4_12_12","P4_12_12",""},/*369*/
{93,"D4^5","P4c2","P4_222","P4_222","P4_222",""},/*370*/
{94,"D4^6","P4n2n","P4_22_12","P4_22_12","P4_22_12",""},/*371*/
{95,"D4^7","P4cw2c","P4_322","P4_322","P4_322",""},/*372*/
{96,"D4^8","P4nw2abw","P4_32_12","P4_32_12","P4_32_12",""},/*373*/
{97,"D4^9","I42","I422","I422","I422",""},/*374*/
{98,"D4^10","I4bw2bw","I4_122","I4_122","I4_122",""},/*375*/
{99,"C4v^1","P4-2","P4mm","P4mm","P4mm",""},/*376*/
{100,"C4v^2","P4-2ab","P4bm","P4bm","P4bm",""},/*377*/
{101,"C4v^3","P4c-2c","P4_2cm","P4_2cm","P4_2cm",""},/*378*/
{102,"C4v^4","P4n-2n","P4_2nm","P4_2nm","P4_2nm",""},/*379*/
{103,"C4v^5","P4-2c","P4cc","P4cc","P4cc",""},/*380*/
{104,"C4v^6","P4-2n","P4nc","P4nc","P4nc",""},/*381*/
{105,"C4v^7","P4c-2","P4_2mc","P4_2mc","P4_2mc",""},/*382*/
{106,"C4v^8","P4c-2ab","P4_2bc","P4_2bc","P4_2bc",""},/*383*/
{107,"C4v^9","I4-2","I4mm","I4mm","I4mm",""},/*384*/
{108,"C4v^10","I4-2c","I4cm","I4cm","I4cm",""},/*385*/
{109,"C4v^11","I4bw-2","I4_1md","I4_1md","I4_1md",""},/*386*/
{110,"C4v^12","I4bw-2c","I4_1cd","I4_1cd","I4_1cd",""},/*387*/
{111,"D2d^1","P-42","P-42m","P-42m","P-42m",""},/*388*/
{112,"D2d^2","P-42c","P-42c","P-42c","P-42c",""},/*389*/
{113,"D2d^3","P-42ab","P-42_1m","P-42_1m","P-42_1m",""},/*390*/
{114,"D2d^4","P-42n","P-42_1c","P-42_1c","P-42_1c",""},/*391*/
{115,"D2d^5","P-4-2","P-4m2","P-4m2","P-4m2",""},/*392*/
{116,"D2d^6","P-4-2c","P-4c2","P-4c2","P-4c2",""},/*393*/
{117,"D2d^7","P-4-2ab","P-4b2","P-4b2","P-4b2",""},/*394*/
{118,"D2d^8","P-4-2n","P-4n2","P-4n2","P-4n2",""},/*395*/
{119,"D2d^9","I-4-2","I-4m2","I-4m2","I-4m2",""},/*396*/
{120,"D2d^10","I-4-2c","I-4c2","I-4c2","I-4c2",""},/*397*/
{121,"D2d^11","I-42","I-42m","I-42m","I-42m",""},/*398*/
{122,"D2d^12","I-42bw","I-42d","I-42d","I-42d",""},/*399*/
{123,"D4h^1","-P42","P4/mmm","P4/m2/m2/m","P4/mmm",""},/*400*/
{124,"D4h^2","-P42c","P4/mcc","P4/m2/c2/c","P4/mcc",""},/*401*/
{125,"D4h^3","P42-1ab","P4/nbm","P4/n2/b2/m","P4/nbm","1"},/*402*/
{125,"D4h^3","-P4a2b","P4/nbm","P4/n2/b2/m","P4/nbm","2"},/*403*/
{126,"D4h^4","P42-1n","P4/nnc","P4/n2/n2/c","P4/nnc","1"},/*404*/
{126,"D4h^4","-P4a2bc","P4/nnc","P4/n2/n2/c","P4/nnc","2"},/*405*/
{127,"D4h^5","-P42ab","P4/mbm","P4/m2_1/bm","P4/mbm",""},/*406*/
{128,"D4h^6","-P42n","P4/mnc","P4/m2_1/nc","P4/mnc",""},/*407*/
{129,"D4h^7","P4ab2ab-1ab","P4/nmm","P4/n2_1/mm","P4/nmm","1"},/*408*/
{129,"D4h^7","-P4a2a","P4/nmm","P4/n2_1/mm","P4/nmm","2"},/*409*/
{130,"D4h^8","P4ab2n-1ab","P4/ncc","P4/n2_1/cc","P4/ncc","1"},/*410*/
{130,"D4h^8","-P4a2ac","P4/ncc","P4/n2_1/cc","P4/ncc","2"},/*411*/
{131,"D4h^9","-P4c2","P4_2/mmc","P4_2/m2/m2/c","P4_2/mmc",""},/*412*/
{132,"D4h^10","-P4c2c","P4_2/mcm","P4_2/m2/c2/m","P4_2/mcm",""},/*413*/
{133,"D4h^11","P4n2c-1n","P4_2/nbc","P4_2/n2/b2/c","P4_2/nbc","1"},/*414*/
{133,"D4h^11","-P4ac2b","P4_2/nbc","P4_2/n2/b2/c","P4_2/nbc","2"},/*415*/
{134,"D4h^12","P4n2-1n","P4_2/nnm","P4_2/n2/n2/m","P4_2/nnm","1"},/*416*/
{134,"D4h^12","-P4ac2bc","P4_2/nnm","P4_2/n2/n2/m","P4_2/nnm","2"},/*417*/
{135,"D4h^13","-P4c2ab","P4_2/mbc","P4_2/m2_1/b2/c","P4_2/mbc",""},/*418*/
{136,"D4h^14","-P4n2n","P4_2/mnm","P4_2/m2_1/n2/m","P4_2/mnm",""},/*419*/
{137,"D4h^15","P4n2n-1n","P4_2/nmc","P4_2/n2_1/m2/c","P4_2/nmc","1"},/*420*/
{137,"D4h^15","-P4ac2a","P4_2/nmc","P4_2/n2_1/m2/c","P4_2/nmc","2"},/*421*/
{138,"D4h^16","P4n2ab-1n","P4_2/ncm","P4_2/n2_1/c2/m","P4_2/ncm","1"},/*422*/
{138,"D4h^16","-P4ac2ac","P4_2/ncm","P4_2/n2_1/c2/m","P4_2/ncm","2"},/*423*/
{139,"D4h^17","-I42","I4/mmm","I4/m2/m2/m","I4/mmm",""},/*424*/
{140,"D4h^18","-I42c","I4/mcm","I4/m2/c2/m","I4/mcm",""},/*425*/
{141,"D4h^19","I4bw2bw-1bw","I4_1/amd","I4_1/a2/m2/d","I4_1/amd","1"},/*426*/
{141,"D4h^19","-I4bd2","I4_1/amd","I4_1/a2/m2/d","I4_1/amd","2"},/*427*/
{142,"D4h^20","I4bw2aw-1bw","I4_1/acd","I4_1/a2/c2/d","I4_1/acd","1"},/*428*/
{142,"D4h^20","-I4bd2c","I4_1/acd","I4_1/a2/c2/d","I4_1/acd","2"},/*429*/
{143,"C3^1","P3","P3","P3","P3",""},/*430*/
{144,"C3^2","P31","P3_1","P3_1","P3_1",""},/*431*/
{145,"C3^3","P32","P3_2","P3_2","P3_2",""},/*432*/
{146,"C3^4","R3","R3","R3","R3","H"},/*433*/
{146,"C3^4","P3*","R3","R3","R3","R"},/*434*/
{147,"C3i^1","-P3","P-3","P-3","P-3",""},/*435*/
{148,"C3i^2","-R3","R-3","R-3","R-3","H"},/*436*/
{148,"C3i^2","-P3*","R-3","R-3","R-3","R"},/*437*/
{149,"D3^1","P32","P312","P312","P312",""},/*438*/
{150,"D3^2","P32=","P321","P321","P321",""},/*439*/
{151,"D3^3","P312c(001)","P3_112","P3_112","P3_112",""},/*440*/
{152,"D3^4","P312=","P3_121","P3_121","P3_121",""},/*441*/
{153,"D3^5","P322c(00-1)","P3_212","P3_212","P3_212",""},/*442*/
{154,"D3^6","P322=","P3_221","P3_221","P3_221",""},/*443*/
{155,"D3^7","R32=","R32","R32","R32","H"},/*444*/
{155,"D3^7","P3*2","R32","R32","R32","R"},/*445*/
{156,"C3v^1","P3-2=","P3m1","P3m1","P3m1",""},/*446*/
{157,"C3v^2","P3-2","P31m","P31m","P31m",""},/*447*/
{158,"C3v^3","P3-2=c","P3c1","P3c1","P3c1",""},/*448*/
{159,"C3v^4","P3-2c","P31c","P31c","P31c",""},/*449*/
{160,"C3v^5","R3-2=","R3m","R3m","R3m","H"},/*450*/
{160,"C3v^5","P3*-2","R3m","R3m","R3m","R"},/*451*/
{161,"C3v^6","R3-2=c","R3c","R3c","R3c","H"},/*452*/
{161,"C3v^6","P3*-2n","R3c","R3c","R3c","R"},/*453*/
{162,"D3d^1","-P32","P-31m","P-312/m","P-31m",""},/*454*/
{163,"D3d^2","-P32c","P-31c","P-312/c","P-31c",""},/*455*/
{164,"D3d^3","-P32=","P-3m1","P-32/m1","P-3m1",""},/*456*/
{165,"D3d^4","-P32=c","P-3c1","P-32/c1","P-3c1",""},/*457*/
{166,"D3d^5","-R32=","R-3m","R-32/m","R-3m","H"},/*458*/
{166,"D3d^5","-P3*2","R-3m","R-32/m","R-3m","R"},/*459*/
{167,"D3d^6","-R32=c","R-3c","R-32/c","R-3c","H"},/*460*/
{167,"D3d^6","-P3*2n","R-3c","R-32/c","R-3c","R"},/*461*/
{168,"C6^1","P6","P6","P6","P6",""},/*462*/
{169,"C6^2","P61","P6_1","P6_1","P6_1",""},/*463*/
{170,"C6^3","P65","P6_5","P6_5","P6_5",""},/*464*/
{171,"C6^4","P62","P6_2","P6_2","P6_2",""},/*465*/
{172,"C6^5","P64","P6_4","P6_4","P6_4",""},/*466*/
{173,"C6^6","P6c","P6_3","P6_3","P6_3",""},/*467*/
{174,"C3h^1","P-6","P-6","P-6","P-6",""},/*468*/
{175,"C6h^1","-P6","P6/m","P6/m","P6/m",""},/*469*/
{176,"C6h^2","-P6c","P6_3/m","P6_3/m","P6_3/m",""},/*470*/
{177,"D6^1","P62","P622","P622","P622",""},/*471*/
{178,"D6^2","P612(00-1)","P6_122","P6_122","P6_122",""},/*472*/
{179,"D6^3","P652(001)","P6_522","P6_522","P6_522",""},/*473*/
{180,"D6^4","P622c(001)","P6_222","P6_222","P6_222",""},/*474*/
{181,"D6^5","P642c(00-1)","P6_422","P6_422","P6_422",""},/*475*/
{182,"D6^6","P6c2c","P6_322","P6_322","P6_322",""},/*476*/
{183,"C6v^1","P6-2","P6mm","P6mm","P6mm",""},/*477*/
{184,"C6v^2","P6-2c","P6cc","P6cc","P6cc",""},/*478*/
{185,"C6v^3","P6c-2","P6_3cm","P6_3cm","P6_3cm",""},/*479*/
{186,"C6v^4","P6c-2c","P6_3mc","P6_3mc","P6_3mc",""},/*480*/
{187,"D3h^1","P-62","P-6m2","P-6m2","P-6m2",""},/*481*/
{188,"D3h^2","P-6c2","P-6c2","P-6c2","P-6c2",""},/*482*/
{189,"D3h^3","P-6-2","P-62m","P-62m","P-62m",""},/*483*/
{190,"D3h^4","P-6c-2c","P-62c","P-62c","P-62c",""},/*484*/
{191,"D6h^1","-P62","P6/mmm","P6/m2/m2/m","P6/mmm",""},/*485*/
{192,"D6h^2","-P62c","P6/mcc","P6/m2/c2/c","P6/mcc",""},/*486*/
{193,"D6h^3","-P6c2","P6_3/mcm","P6_3/m2/c2/m","P6_3/mcm",""},/*487*/
{194,"D6h^4","-P6c2c","P6_3/mmc","P6_3/m2/m2/c","P6_3/mmc",""},/*488*/
{195,"T^1","P223","P23","P23","P23",""},/*489*/
{196,"T^2","F223","F23","F23","F23",""},/*490*/
{197,"T^3","I223","I23","I23","I23",""},/*491*/
{198,"T^4","P2ac2ab3","P2_13","P2_13","P2_13",""},/*492*/
{199,"T^5","I2b2c3","I2_13","I2_13","I2_13",""},/*493*/
{200,"Th^1","-P223","Pm3","P2/m-3","Pm3",""},/*494*/
{201,"Th^2","P223-1n","Pn3","P2/n-3","Pn3","1"},/*495*/
{201,"Th^2","-P2ab2bc3","Pn3","P2/n-3","Pn3","2"},/*496*/
{202,"Th^3","-F223","Fm3","F2/m-3","Fm3",""},/*497*/
{203,"Th^4","F223-1d","Fd3","F2/d-3","Fd3","1"},/*498*/
{203,"Th^4","-F2uv2vw3","Fd3","F2/d-3","Fd3","2"},/*499*/
{204,"Th^5","-I223","Im3","I2/m-3","Im3",""},/*500*/
{205,"Th^6","-P2ac2ab3","Pa3","P2_1/a-3","Pa3",""},/*501*/
{206,"Th^7","-I2b2c3","Ia3","I2_1/a-3","Ia3",""},/*502*/
{207,"O^1","P423","P432","P432","P432",""},/*503*/
{208,"O^2","P4n23","P4_232","P4_232","P4_232",""},/*504*/
{209,"O^3","F423","F432","F432","F432",""},/*505*/
{210,"O^4","F4d23","F4_132","F4_132","F4_132",""},/*506*/
{211,"O^5","I423","I432","I432","I432",""},/*507*/
{212,"O^6","P4acd2ab3","P4_332","P4_332","P4_332",""},/*508*/
{213,"O^7","P4bd2ab3","P4_132","P4_132","P4_132",""},/*509*/
{214,"O^8","I4bd2c3","I4_132","I4_132","I4_132",""},/*510*/
{215,"Td^1","P-423","P-43m","P-43m","P-43m",""},/*511*/
{216,"Td^2","F-423","F-43m","F-43m","F-43m",""},/*512*/
{217,"Td^3","I-423","I-43m","I-43m","I-43m",""},/*513*/
{218,"Td^4","P-4n23","P-43n","P-43n","P-43n",""},/*514*/
{219,"Td^5","F-4c23","F-43c","F-43c","F-43c",""},/*515*/
{220,"Td^6","I-4bd2c3","I-43d","I-43d","I-43d",""},/*516*/
{221,"Oh^1","-P423","Pm-3m","P4/m-32/m","Pm-3m",""},/*517*/
{222,"Oh^2","P423-1n","Pn-3n","P4/n-32/n","Pn-3n","1"},/*518*/
{222,"Oh^2","-P4a2bc3","Pn-3n","P4/n-32/n","Pn-3n","2"},/*519*/
{223,"Oh^3","-P4n23","Pm-3n","P4_2/m-32/n","Pm-3n",""},/*520*/
{224,"Oh^4","P4n23-1n","Pn-3m","P4_2/n-32/m","Pn-3m","1"},/*521*/
{224,"Oh^4","-P4bc2bc3","Pn-3m","P4_2/n-32/m","Pn-3m","2"},/*522*/
{225,"Oh^5","-F423","Fm-3m","F4/m-32/m","Fm-3m",""},/*523*/
{226,"Oh^6","-F4c23","Fm-3c","F4/m-32/c","Fm-3c",""},/*524*/
{227,"Oh^7","F4d23-1d","Fd-3m","F4_1/d-32/m","Fd-3m","1"},/*525*/
{227,"Oh^7","-F4vw2vw3","Fd-3m","F4_1/d-32/m","Fd-3m","2"},/*526*/
{228,"Oh^8","F4d23-1cd","Fd-3c","F4_1/d-32/c","Fd-3c","1"},/*527*/
{228,"Oh^8","-F4cvw2vw3","Fd-3c","F4_1/d-32/c","Fd-3c","2"},/*528*/
{229,"Oh^9","-I423","Im-3m","I4/m-32/m","Im-3m",""},/*529*/
{230,"Oh^10","-I4bd2c3","Ia-3d","I4_1/a-32/d","Ia-3d",""},/*530*/
};

// Process known handles for all windows
void guiapps_process_handles(){
  // close all windows?
  if (settings.close_all_windows){
    structureinfo_window_h = false;
    structurenew_window_h = false;
    structureopen_window_h = false;
    console_window_h = false;
    if (settings.preview_mode)
      structurenew_window(&structurenew_window_h);
    settings.close_all_windows = false;
  }

  // window handles
  if (structureinfo_window_h) structureinfo_window(&structureinfo_window_h);
  if (structurenew_window_h) structurenew_window(&structurenew_window_h);
  if (structureopen_window_h > 0) structureopen_window(&structureopen_window_h);
  if (console_window_h) console_window(&console_window_h);

  // overlay showing the current mode
  if (settings.preview_mode){
    ImGui::SetNextWindowPos(ImVec2(10,30));
    if (ImGui::Begin("###modeoverlay", NULL, ImVec2(0,0), 0.3f, ImGuiWindowFlags_NoTitleBar|
		     ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoSavedSettings)){
      ImGui::TextColored(ImVec4(1.0f,1.0f,0.0f,1.0f),"Preview");
    }
    ImGui::End();
  }
}

// Open a layout showing the structural information for this 
// crystal/molecule
void structureinfo_window(bool *p_open){
  const char *title[] = {"Structure not available###structinfo", "Molecule information###structinfo", "Crystal information###structinfo"};

  const int elem_unselected = 0;
  const int elem_general_unavailable = 1;
  const int elem_general_molecule = 2;
  const int elem_general_crystal = 3;
  const int asym_unit_crystal = 4;
  const int unit_cell_crystal = 5;
  const int symmetry_crystal = 6;
  const int fragments_crystal = 7;
  const int atoms_molecule = 8;
  const int fragments_molecule = 9;

  const int nelems = 10;
  const int elem_levels[nelems] = {0,0,1,2,2,2,2,2,1,1};

  const char *elem_title[] = {"", "", "General","General","Asymmetric unit","Unit cell",
			      "Symmetry","Fragments","Atoms","Fragments"};

  // Display different information depending on the type of structure (level)
  int level;
  if (!isinit)
    level = 0;
  else if (ismolecule)
    level = 1;
  else
    level = 2;

  static int selected = 0;

  ImGui::SetNextWindowSize(ImVec2(640, 350), ImGuiSetCond_Once);
  ImGui::SetNextWindowPos(ImVec2(0, 20), ImGuiSetCond_Once);
  if (ImGui::Begin(title[level], p_open)){
    ImGui::BeginGroup();

    // List of selectable elements on the left
    ImGui::BeginChild("left", ImVec2(150, 0), true);
    for (int i=2;i<nelems;i++){
      if (elem_levels[i] == level && ImGui::Selectable(elem_title[i], selected == i)){
	selected = i;
      }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    if (elem_levels[selected] != level)
      selected = 0;

    ImGui::BeginChild("right", ImVec2(0, -ImGui::GetItemsLineHeightWithSpacing())); 
    ImGui::Text(elem_title[selected]);
    ImGui::Separator();

    char *text = get_text_info(selected);
    if (text){
      ImGui::TextWrapped(text);
      free(text);
    }
    ImGui::EndChild();

    ImGui::EndGroup();
  }
  ImGui::End();
}

// Create a new structure with user's input
void structurenew_window(bool *p_open){
  static struct c_crystalseed useed = {-1};
  static Settings *settings0;

  // xxxx error incorrect space group if empty get current and then molecule
  // xxxx use the spg chooser
  // xxxx mgo 5 bonds?
  // xxxx delete and esc keybindings
  // xxxx cartesian axes
  // xxxx molecular cell

  static int atunitsc, atunitsm;
  const char *latstr[] = {"Lengths & angles","Lattice vectors"};
  const char *unitstr[] = {"Bohr","Angstrom","Fractional"};
  bool doclear = false;

  if (useed.type == -1){
    useed.molborder = 10.0f * 0.529177f;
    useed.borunits = 1;
    atunitsc = 2;
    atunitsm = 1;
  }

  ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiSetCond_Once);
  ImGui::SetNextWindowSizeConstraints(ImVec2(500, 450), ImVec2(FLT_MAX, FLT_MAX)); 
  ImGui::SetNextWindowPos(ImVec2(150, 150), ImGuiSetCond_Once);
  if (ImGui::Begin("New structure", p_open)){
    ImGui::BeginChild("###appbody", ImVec2(0, -2*ImGui::GetItemsLineHeightWithSpacing()));
    // structure type
    ImGui::RadioButton("Crystal", &useed.type, 1);
    ImGui::SameLine();
    ImGui::RadioButton("Molecule", &useed.type, 0);
    ImGui::SameLine();
    if (ImGui::Button("Get Current")){
      useed = get_seed_from_current_structure();
      if (settings0 && useed.type != -1){
	reject_previewed_structure();
	settings = *settings0;
	delete settings0;
	settings0 = NULL;
	settings.preview_mode = false;
      }
    }
    ImGui::SameLine();
    doclear = ImGui::Button("Clear");

    if (useed.type == 1){
      ImGui::Separator();
      ImGui::Text("Cell");
      ImGui::SameLine();
      ImGui::PushItemWidth(170);
      ImGui::Combo("###achoice", &useed.achoice, latstr, 2);
      ImGui::PopItemWidth();

      if (useed.achoice == 0){
	// Cell lengths and angles
	ImGui::Text("Lengths");
	ImGui::SameLine();
	ImGui::InputText("###straa",useed.straa,IM_ARRAYSIZE(useed.straa));
	ImGui::SameLine();
	ImGui::PushItemWidth(90);
	ImGui::Combo("###aaunits", &useed.aaunits, unitstr, 2);
	ImGui::PopItemWidth();

	ImGui::Text("Angles ");
	ImGui::SameLine();
	ImGui::InputText("###strbb",useed.strbb,IM_ARRAYSIZE(useed.strbb));
	ImGui::SameLine();
	ImGui::Text("degrees");
      }
      else{
	// Cartesian coordinates of the lattice vectors
	ImGui::Text("Lattice vectors:");
	ImGui::SameLine();
	ImGui::PushItemWidth(90);
	ImGui::Combo("###rrunits", &useed.rrunits, unitstr, 2);
	ImGui::PopItemWidth();
	ImGui::InputTextMultiline("###strrr", useed.strrr, IM_ARRAYSIZE(useed.strrr), ImVec2(-1.0f, ImGui::GetTextLineHeight() * 5));
	ImGui::Text("(each vector in a row; arithmetic expressions allowed)");
      }

      // Space group
      ImGui::Separator();
      ImGui::Text("Space group: ");
      ImGui::SameLine();
      ImGui::PushItemWidth(120);
      ImGui::InputText("###strspg",useed.strspg,IM_ARRAYSIZE(useed.strspg));
      ImGui::PopItemWidth();
      ImGui::SameLine();
      if (ImGui::Button("Choose it")){}
      // xxxx spg chooser xxxx //
      // if (ImGui::Button("Choose it"))
      // 	ImGui::OpenPopup("choosespg");
      // if (ImGui::BeginPopup("choosespg")){
      // 	int lspgnum = 0;
      // 	if (ImGui::BeginMenu("Number")){
      // 	  lspgnum = spg_choose_menu(0);
      // 	  ImGui::EndMenu();
      // 	}
      // 	if (ImGui::BeginMenu("Hermann-Mauguin")){
      // 	  lspgnum = spg_choose_menu(1);
      // 	  ImGui::EndMenu();
      // 	}
      // 	if (ImGui::BeginMenu("Hall")){
      // 	  lspgnum = spg_choose_menu(2);
      // 	  ImGui::EndMenu();
      // 	}
      // 	if (ImGui::BeginMenu("Schoenflies")){
      // 	  lspgnum = spg_choose_menu(3);
      // 	  ImGui::EndMenu();
      // 	}
      // 	if (lspgnum > 0){
      // 	  strcpy(useed.strspg,spgtyp[lspgnum].international_full);
      // 	}
      // 	ImGui::EndPopup();
      // }
    }
    
    if (useed.type > -1){
      // Atomic positions
      ImGui::Separator();
      ImGui::Text("Atomic coordinates:");
      ImGui::SameLine();
      ImGui::PushItemWidth(110);
      if (useed.type == 1){
	if (useed.achoice == 0){
	  ImGui::Text(" (Fractional)");
	  useed.atunits = 2;
	}
	else{
	  ImGui::Combo("###atunits", &atunitsc, unitstr, 3);
	  useed.atunits = atunitsc;
	}
      }
      else {
	ImGui::Combo("###atunits", &atunitsm, unitstr, 2);
	useed.atunits = atunitsm;
      }
      ImGui::PopItemWidth();
      ImGui::InputTextMultiline("###strat", useed.strat, IM_ARRAYSIZE(useed.strat),
				ImVec2(-1.0f, ImGui::GetTextLineHeight() * 10));
      ImGui::Text("(each atom in a row: symbol x y z; arithmetic expressions allowed)");
    }

    if (useed.type == 0){
      ImGui::Checkbox("Use a cubic cell", &useed.molcubic);

      ImGui::Text("Cell border:");
      ImGui::SameLine();
      ImGui::PushItemWidth(90);
      ImGui::DragFloat("###bordrag", &useed.molborder, 0.02f, 0.f, FLT_MAX, "%.3f");
      ImGui::SameLine();
      ImGui::Combo("###borunits", &useed.borunits, unitstr, 2);
      ImGui::PopItemWidth();
    }
    ImGui::EndChild();

    // Error message at the end of the window but before the buttons
    ImGui::BeginChild("###errmsg", ImVec2(0, ImGui::GetItemsLineHeightWithSpacing()));
    if (useed.errcode)
      ImGui::TextColored(ImVec4(1.0f,0.0f,0.0f,1.0f), "Error: %s",useed.errmsg);
    ImGui::EndChild();

    // Buttons at the bottom of the window
    ImGui::BeginChild("###buttons", ImVec2(0, ImGui::GetItemsLineHeightWithSpacing()));
    ImGui::PushItemWidth(-140);
    ImGui::SameLine(ImGui::GetWindowContentRegionWidth()-150);
    if (ImGui::Button("Preview")){
      if (preview_structure(&useed)){
	// save a copy of the old settings
	if (!settings0) {
	  settings0 = new Settings(ismolecule,box_xmaxlen);
	  *settings0 = settings;
	}
	// settings for the previewed crystal
	settings.set_flags_and_cam(useed.type == 0,box_xmaxlen,box_xmaxclen);
	settings.preview_mode = true;
      } 
    }

    ImGui::SameLine();
    if (ImGui::Button("OK")) {
      if (new_structure(&useed,false)){
	if (settings0){
	  accept_previewed_structure();
	  delete settings0;
	  settings0 = NULL;
	}
	else{
	  settings.set_flags_and_cam(useed.type == 0,box_xmaxlen,box_xmaxclen);
	}
	*p_open = false;
	settings.preview_mode = false;
      } 
    }

    ImGui::SameLine();
    if (ImGui::Button("Cancel")) 
      *p_open = false;
    ImGui::PopItemWidth();
    ImGui::EndChild();
  }
  ImGui::End();

  // reset the static variables for the next use
  if (!*p_open || doclear){
    useed = {-1};
    if (settings0){
      reject_previewed_structure();
      settings = *settings0;
      delete settings0;
      settings0 = NULL;
      settings.preview_mode = false;
    }
  }
}

// Space group choosing menu ietms (beginmenu and endmenu must be
// provided). mode = 0 (number), 1 (HM), 2 (Hall), 3 (Schoenflies).
static int spg_choose_menu(int mode){
  bool spgflag;
  int n = 1;
  char label[128];

  ImGui::Columns(5,"###spgnumber",true);
  for (int i=0;i<5;i++){
    for (int j=0;j<106;j++,n++){
      switch (mode){
      case 0:
	if (strlen(spgtyp[n].choice) > 0)
	  sprintf(label, "%d (%s)",spgtyp[n].number,spgtyp[n].choice);
	else
	  sprintf(label, "%d",spgtyp[n].number);
	break;
      case 1:
	if (strlen(spgtyp[n].choice) > 0)
	  sprintf(label, "%s (%s)",spgtyp[n].international_full,spgtyp[n].choice);
	else
	  sprintf(label, "%s",spgtyp[n].international_full);
	break;
      case 2:
	sprintf(label, "%s",spgtyp[n].hall_symbol);
	break;
      case 3:
	if (strlen(spgtyp[n].choice) > 0)
	  sprintf(label, "%s (%s)",spgtyp[n].schoenflies,spgtyp[n].choice);
	else
	  sprintf(label, "%s",spgtyp[n].schoenflies);
	break;
      }
      if (ImGui::MenuItem(label,NULL,false))
	return n;
    }
    ImGui::NextColumn();
  }
  ImGui::Columns(1);
  return 0;
}

// Read the structure from an external file p_open = 0 - closed, 1 - molecule, 2 - crystal
void structureopen_window(int *p_open){
  static ImGuiFs::Dialog fsopenfile;
  static bool firstpass = true;

  const char* filename = fsopenfile.chooseFileDialog(firstpass,"./",NULL);
  firstpass = false;

  if (fsopenfile.hasUserJustCancelledDialog() && strlen(filename) == 0){
    // Dialog has been closed - set up for next time and prevent more calls for now
    firstpass = true;
    *p_open = 0;
  }

  if (strlen(filename) > 0){
    // Clean up previous and initialize the structure
    open_structure(&filename, *p_open == 1); 
    settings.set_flags_and_cam(*p_open == 1,box_xmaxlen,box_xmaxclen);

    // Close the dialog
    firstpass = true;
    *p_open = 0;
  }
}

void console_window(bool *p_open){
  static char command[2048] = "";

  ImGuiIO& io = ImGui::GetIO();

  ImGui::SetNextWindowPos(ImVec2(10,io.DisplaySize.y-2*ImGui::GetTextLineHeightWithSpacing()));
  ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x,ImGui::GetTextLineHeightWithSpacing()));
  if (ImGui::Begin("",NULL,ImVec2(0.f,0.f),0.0,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoSavedSettings)){
    ImGui::Text("Input:");
    ImGui::SameLine();
    ImGui::SetKeyboardFocusHere();
    if (ImGui::InputText("###inputconsole", command, IM_ARRAYSIZE(command), ImGuiInputTextFlags_EnterReturnsTrue|ImGuiInputTextFlags_AutoSelectAll|ImGuiInputTextFlags_AlwaysInsertMode)){
      run_critic2_command(command);
      command[0] = '\0';
      *p_open = false;
    }
    ImGui::End();
  }
}
