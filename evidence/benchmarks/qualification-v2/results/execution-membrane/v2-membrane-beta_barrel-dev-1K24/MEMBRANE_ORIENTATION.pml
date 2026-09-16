# memorient oriented view
# context=gram_negative_om method=barrel_normal label=surface
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 5+6+7+8+43+44+45+46+47+48+100+101+102+103+104+105+151+202+203+204
color acc_periplasmic, periplasmic
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 9+11+13+15+39+41+51+53+91+93+97+99+107+109+138+140+142+144+146+148+154+156+158+160+162+196+200+206+208+212+248+252
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 10+12+14+36+37+38+40+42+49+50+52+54+55+56+90+92+94+95+96+98+106+108+110+111+112+113+114+115+116+137+139+141+143+145+147+149+150+152+153+155+157+159+161+163+165+191+192+193+194+195+197+198+199+201+205+207+209+210+211+247+249+250+251+253
color acc_lipid_embedded, lipid_embedded
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 16+17+18+19+32+33+34+35+57+58+59+86+87+88+89+117+118+119+120+133+134+135+136+164+166+167+168+185+186+187+188+189+190+213+214+215+216+217+218+242+243+244+245+246
color acc_lps_shielded, lps_shielded
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 20+21+30+31+66+68+70+71+72+82+84+85+121+123+125+128+130+132+169+170+171+172+173+176+178+182+183+184+219+221+224+226+228+233+235+237+238+240
color acc_buried_interior, buried_interior
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 22+23+24+25+26+27+28+29+60+61+62+63+64+65+67+69+73+74+75+76+77+78+79+80+81+83+122+124+126+127+129+131+174+175+177+179+180+181+220+222+223+225+227+229+230+231+232+234+236+239+241
color acc_antibody_accessible, antibody_accessible
select epitope_surface, mol and resi 22+23+24+25+26+27+28+29+60+61+62+63+64+65+67+69+73+74+75+76+77+78+79+80+81+83+122+124+126+127+129+131+174+175+177+179+180+181+220+222+223+225+227+229+230+231+232+234+236+239+241
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,11.1]
pseudoatom mem_peri, pos=[0,0,-11.1]
# membrane core spans z = [-11.1, 11.1]; extracellular is +z
orient mol
zoom mol, 5
