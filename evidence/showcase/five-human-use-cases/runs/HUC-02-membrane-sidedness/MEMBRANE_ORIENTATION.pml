# memorient oriented view
# context=tm_receptor method=tm_helix_axis_v2 label=tm_helix_experimental
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 1+2+3+4+5+6+7
color acc_periplasmic, periplasmic
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29
color acc_lipid_embedded, lipid_embedded
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57
color acc_antibody_accessible, antibody_accessible
select epitope_surface, mol and resi 30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,16.0]
pseudoatom mem_peri, pos=[0,0,-16.0]
# membrane core spans z = [-16.0, 16.0]; extracellular is +z
orient mol
zoom mol, 5
