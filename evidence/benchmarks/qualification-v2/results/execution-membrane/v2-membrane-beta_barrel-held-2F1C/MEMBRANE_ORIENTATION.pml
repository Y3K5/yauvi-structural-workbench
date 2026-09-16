# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 4+5+7+9+13+35+37+39+44+46+48+72+74+76+79+84+88+115+117+119+121+127+129+131+133+152+154+158+160+166+168+170+172+192+194+196+198+200+205+207+209+211+235+239+243+245+249+251+253+255+277+279+281
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 6+8+10+11+12+34+36+38+40+41+42+43+45+47+71+73+75+77+78+83+85+86+87+89+90+114+116+118+120+122+126+128+130+132+153+155+156+157+159+165+167+169+171+173+193+195+197+199+201+203+204+206+208+210+212+213+214+234+236+237+238+240+241+242+244+246+247+248+250+252+254+256+272+273+274+275+276+278+280+282+283
color acc_lipid_embedded, lipid_embedded
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 14+15+16+17+18+19+29+30+31+32+33+49+50+51+52+67+68+69+70+91+92+93+110+111+112+113+134+135+136+149+150+151+174+175+176+188+189+190+191+215+216+217+232+233+257+258+259+260+267+268+269+270+271
color acc_lps_shielded, lps_shielded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 28+55+56+57+58+61+62+64+100+101+102+103+104+105+106+108+139+142+179+180+181+182+184+218+219
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 53+54+63+65+66+94+95+96+97+98+99+107+109+137+138+140+141+143+144+145+146+147+148+177+178+183+185+186+187
color acc_buried_interior, buried_interior
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 80+81+82+123+124+125+161+162+163+164+202
color acc_periplasmic, periplasmic
select epitope_surface, mol and resi 28+55+56+57+58+61+62+64+100+101+102+103+104+105+106+108+139+142+179+180+181+182+184+218+219
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.0]
pseudoatom mem_peri, pos=[0,0,-12.0]
# membrane core spans z = [-12.0, 12.0]; extracellular is +z
orient mol
zoom mol, 5
