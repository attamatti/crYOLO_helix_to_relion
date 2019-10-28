Convert the EMAN format boxes from crYOLO to relion format while preserving the helical tube IDs for angular priors

Instructions:
1) run crYOLO_helix_parts.star 
	USAGE: crYOLO_helix_parts.star <crYOLO boxfile dir>

2) Use relion to extract the segments as single particles (not helical segments)
	use the extract function with coords_suffix_YOLObox.box (which is 1 dir up from the crYOLO bosfile dir) as the particles coordinates 

3) run rln3p1_crYOLO_add_filaments.py to add the helical tube IDs to the particles file
	USAGE: rln3p1_crYOLO_add_filaments.py <extracted paricles star file> <crYOLO boxfile dir>

4) Buy Matt some cake and give an acknowledgement in you paper 
