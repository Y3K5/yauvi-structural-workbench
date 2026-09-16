# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 1+3+5+7+19+20+21+23+25+37+39+41+43+44+45+60+62+64+66+68+70+75+77+79+81+82+83+84+85+104+106+108+109+110+111+112+114+115+116+117+118+120+121+122+124+142+144+146+148+149+150
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 2+4+6+18+22+24+36+38+40+42+46+61+63+65+67+69+76+78+80+86+105+107+113+119+123+125+143+145+147
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 8+9+10+11+13+14+15+16+17+47+48+49+56+57+58+59+87+88+89+90+98+100+101+102+103+126+127+128+129+130+131+132+133+136+138+139+140+141
color acc_lps_shielded, lps_shielded
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 12+50+54+91+93+94+95+96+99+135+137
color acc_buried_interior, buried_interior
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 26+27+28+29+30+31+32+33+34+35+71+72+73+74
color acc_periplasmic, periplasmic
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 51+52+53+55+92+97+134
color acc_antibody_accessible, antibody_accessible
select epitope_surface, mol and resi 51+52+53+55+92+97+134
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.6]
pseudoatom mem_peri, pos=[0,0,-12.6]
# membrane core spans z = [-12.6, 12.6]; extracellular is +z
orient mol
zoom mol, 5
