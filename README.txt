Convert the EMAN format boxes from crYOLO to relion format while preserving the helical tube IDs for angular priors

Instructions:
1) run rln_crYOLO_helix_box.py 
	USAGE: rln_crYOLO_helix_box.star <crYOLO boxfile dir>
	Note the crYOLO boxfile dir needs to have the same name as the dir with your raw data

2) Use relion to extract the segments as single particles (not helical segments)
	use the extract function with coords_suffix_YOLObox.box (which is 1 dir up from the crYOLO bosfile dir) as the particles' coordinates 

3) run rln3p1_crYOLO_helix_write_priors.py to add the helical tube IDs, psi priors, and helical track lengths to the particles file
	USAGE: rln3p1_crYOLO_helix_write_priors.py <extracted paricles star file> <crYOLO boxfile dir> <amout of overalp Angstrom> <px size Angstrom> <box size pixels>

4) Buy Matt some cake and give an acknowledgement in your paper 
