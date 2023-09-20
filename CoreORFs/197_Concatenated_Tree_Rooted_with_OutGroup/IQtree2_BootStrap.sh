#!/bin/bash

iqtree2 \
	-T AUTO \
	-B 1000 \
	-mset raxml \
	--threads-max 50 \
	-o Clostridium_thermobutyricum \
	-s 197_OrthoGroups_Alignment.fa \
	-pre 197_OrthoGroups_Alignment \
	

# iqtree2 \
	# -T AUTO \
	# -B 1000 \
	# -mset raxml \
	# --threads-max 50 \
	# -o Clostridium_thermobutyricum \
	# -s 01.Alignment/ \
	# -pre 197_OrthoGroups_Alignment \
	

