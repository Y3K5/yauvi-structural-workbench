# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 9+10+11+12+45+46+88+91+92+93+134+136+175+178+225+273+274+275
color acc_periplasmic, periplasmic
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 13+15+17+19+21+23+35+37+39+41+43+47+48+50+52+54+56+57+58+76+78+80+82+84+86+90+94+96+97+98+99+102+104+106+108+110+121+123+125+126+127+128+129+131+132+135+137+139+141+143+145+147+163+165+167+169+171+173+174+176+177+179+181+183+185+187+189+214+216+218+220+222+224+226+227+229+231+233+235+262+263+264+266+268+269+270+272
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 14+16+18+20+22+34+36+38+40+42+44+49+51+53+55+75+77+79+81+83+85+87+89+95+100+101+103+105+107+109+111+122+124+130+133+138+140+142+144+146+164+166+168+170+172+180+182+184+186+188+190+215+217+219+221+223+228+230+232+234+236+265+267+271
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 24+25+26+27+28+29+30+31+32+33+112+113+114+115+116+117+118+119+120+148+149+150+151+152+156+157+158+159+160+161+162+191+192+193+194+195+210+211+212+213+237+238+239+240+241+242+244+258+259+260+261
color acc_lps_shielded, lps_shielded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 153+155+197+198+200+201+203+204+205+206+207+245+249+250+251+252+253+254+256
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 154+196+199+202+208+209+243+246+247+248+255+257
color acc_buried_interior, buried_interior
select epitope_surface, mol and resi 153+155+197+198+200+201+203+204+205+206+207+245+249+250+251+252+253+254+256
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,13.5]
pseudoatom mem_peri, pos=[0,0,-13.5]
# membrane core spans z = [-13.5, 13.5]; extracellular is +z
orient mol
zoom mol, 5
