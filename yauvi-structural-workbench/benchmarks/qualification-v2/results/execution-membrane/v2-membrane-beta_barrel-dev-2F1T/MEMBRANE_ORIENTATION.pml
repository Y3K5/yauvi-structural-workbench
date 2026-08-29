# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 2+3+48+49+90+91+92+93+94+95+145+146
color acc_periplasmic, periplasmic
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 4+6+8+10+11+12+14+39+40+41+43+45+47+50+51+52+53+55+57+81+83+85+87+89+96+98+100+101+102+103+104+133+134+136+138+139+140+142+144+147+148+150+152+154+156+181+183+185+187+188+189+191
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 5+7+9+13+42+44+46+54+56+58+82+84+86+88+97+99+105+135+137+141+143+149+151+153+155+157+182+184+186+190+192
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 15+16+17+18+19+20+34+35+36+37+38+59+60+61+62+63+77+78+79+80+106+107+108+109+130+131+132+158+159+160+161+177+178+179+180
color acc_lps_shielded, lps_shielded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 29+30+32+64+68+69+70+71+72+75+110+111+113+115+116+117+119+120+121+122+124+125+127+129+162+166+167+168+169+170+171+173+174+175
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 31+33+65+66+67+73+74+76+112+114+118+123+126+128+163+164+165+172+176
color acc_buried_interior, buried_interior
select epitope_surface, mol and resi 29+30+32+64+68+69+70+71+72+75+110+111+113+115+116+117+119+120+121+122+124+125+127+129+162+166+167+168+169+170+171+173+174+175
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.5]
pseudoatom mem_peri, pos=[0,0,-12.5]
# membrane core spans z = [-12.5, 12.5]; extracellular is +z
orient mol
zoom mol, 5
