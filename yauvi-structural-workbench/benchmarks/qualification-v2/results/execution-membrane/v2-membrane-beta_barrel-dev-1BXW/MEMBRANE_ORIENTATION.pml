# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 0+2+3+4+5+6+46+47+88+89+90+132+133+134+170+171
color acc_periplasmic, periplasmic
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 1+7+9+10+11+13+14+37+38+39+41+42+43+45+48+49+50+51+53+54+55+76+77+79+81+83+84+85+87+91+93+95+97+98+99+101+121+123+125+126+127+129+131+135+137+139+141+143+162+164+165+166+168
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 8+12+40+44+52+56+78+80+82+86+92+94+96+100+122+124+128+130+136+138+140+142+144+161+163+167+169
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 15+16+17+18+19+33+34+35+36+57+58+59+60+61+73+74+75+102+103+104+117+118+119+120+145+146+147+148+149+156+157+158+159+160
color acc_lps_shielded, lps_shielded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 20+21+22+23+24+25+26+27+28+29+30+31+62+63+64+65+66+67+68+69+70+71+105+106+107+108+109+110+111+112+113+114+115+116+150+151+152+153+154+155
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 32+72
color acc_buried_interior, buried_interior
select epitope_surface, mol and resi 20+21+22+23+24+25+26+27+28+29+30+31+62+63+64+65+66+67+68+69+70+71+105+106+107+108+109+110+111+112+113+114+115+116+150+151+152+153+154+155
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.3]
pseudoatom mem_peri, pos=[0,0,-12.3]
# membrane core spans z = [-12.3, 12.3]; extracellular is +z
orient mol
zoom mol, 5
