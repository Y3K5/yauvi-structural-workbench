# memorient oriented view
# context=gram_negative_om method=barrel_normal label=barrel
load oriented.pdb, mol
hide everything, mol
show cartoon, mol
bg_color white
set cartoon_transparency, 0.1
set_color acc_antibody_accessible, [0.902, 0.624, 0.000]
select antibody_accessible, mol and resi 786+788+833+835+836+838+840+868+870+872+874+875+877+879+881+907+909+913+915+916+917+918+919+920+921+923+925+961+963+964+966+968+969+971+972+973+975+977+978+980+1016+1018+1020+1021
color acc_antibody_accessible, antibody_accessible
set_color acc_buried_interior, [0.600, 0.600, 0.600]
select buried_interior, mol and resi 787+789+790+791+792+793+794+795+796+797+834+837+839+841+869+871+873+876+878+880+882+908+910+911+912+914+922+924+926+927+928+929+960+962+965+967+970+974+976+979+981+1015+1017+1019
color acc_buried_interior, buried_interior
set_color acc_lps_shielded, [0.337, 0.706, 0.914]
select lps_shielded, mol and resi 798+799+800+801+802+803+804+826+827+828+829+830+831+832+842+843+844+845+846+847+864+865+866+867+883+884+885+886+904+905+906+930+931+932+933+957+958+959+982+983+984+985+1010+1011+1012+1013+1014+1039+1040+1041+1067+1068+1069+1070+1071+1072+1073+1074
color acc_lps_shielded, lps_shielded
set_color acc_lipid_embedded, [0.941, 0.894, 0.259]
select lipid_embedded, mol and resi 805+807+808+811+812+815+816+817+818+819+821+823+825+848+849+850+852+854+855+856+857+858+860+862+887+888+889+891+893+894+896+897+899+901+903+934+935+937+938+939+941+947+949+951+953+954+955+986+988+989+990+992+994+996+999+1000+1002+1004+1006+1007+1008+1043+1045+1046+1047+1048+1049+1051+1053+1054+1055+1056+1057+1059+1061+1063+1065+1076+1078+1079+1080+1081+1082+1084
color acc_lipid_embedded, lipid_embedded
set_color acc_pore_lumen_facing, [0.000, 0.447, 0.698]
select pore_lumen_facing, mol and resi 806+809+810+813+814+820+822+824+851+853+859+861+863+890+892+895+898+900+902+936+940+948+950+952+956+987+991+993+995+1001+1003+1005+1009+1042+1044+1050+1052+1058+1060+1062+1064+1066+1075+1077+1083
color acc_pore_lumen_facing, pore_lumen_facing
set_color acc_periplasmic, [0.835, 0.369, 0.000]
select periplasmic, mol and resi 942+946+997+998
color acc_periplasmic, periplasmic
select epitope_surface, mol and resi 786+788+833+835+836+838+840+868+870+872+874+875+877+879+881+907+909+913+915+916+917+918+919+920+921+923+925+961+963+964+966+968+969+971+972+973+975+977+978+980+1016+1018+1020+1021
show sticks, epitope_surface
pseudoatom mem_ec, pos=[0,0,12.6]
pseudoatom mem_peri, pos=[0,0,-12.6]
# membrane core spans z = [-12.6, 12.6]; extracellular is +z
orient mol
zoom mol, 5
