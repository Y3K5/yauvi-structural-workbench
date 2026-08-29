# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 1+2+4+38+39+73+74+75+76+77+119+120+121
color acc_periplasmic, periplasmic
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 3+5+6+8+10+12+14+29+31+33+34+35+37+40+42+44+46+64+65+66+68+70+72+78+80+82+83+84+86+88+107+108+109+110+112+113+114+116+118+122+124+126+127+128+130+145+147+149+150+151+153+155
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 7+9+11+13+15+28+30+32+36+41+43+45+47+63+67+69+71+79+81+85+87+111+115+117+123+125+129+131+146+148+152+154
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 16+17+18+19+23+24+25+26+27+48+49+50+51+52+60+61+62+89+90+91+92+103+104+105+106+132+133+134+135+136+141+142+143+144
color acc_lps_shielded, lps_shielded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 20+21+53+54+55+56+57+58+94+95+96+97+98+99+100+101+137+138+140
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 22+59+93+102+139
color acc_buried_interior, buried_interior
select epitope_surface, mol and resi 20+21+53+54+55+56+57+58+94+95+96+97+98+99+100+101+137+138+140
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.7]
pseudoatom mem_peri, pos=[0,0,-12.7]
# membrane core spans z = [-12.7, 12.7]; extracellular is +z
orient mol
zoom mol, 5
