# ---------------------------------------------------------------------
# Data for Fig 1/2 was made in their respective notebooks 
# (in the figs/ dir)
#
# ---------------------------------------------------------------------
fig3: fig3a fig3b fig3c

fig3a:
	-mkdir data/fig3
	nice -19 python exp/fig3.py data/fig3 pars/fig3 a

fig3b:
	-mkdir data/fig3
	nice -19 python exp/fig3.py data/fig3 pars/fig3 b

fig3c:
	-mkdir data/fig3
	nice -19 python exp/fig3.py data/fig3 pars/fig3 c

# ---------------------------------------------------------------------
fig4_lowthres: fig4a_part1 fig4b_part1 fig4a_part4 fig4b_part4


fig4_medthres: fig4a_part2 fig4b_part2 fig4a_part5 fig4b_part5


fig4_highthres: fig4a_part3 fig4b_part3 fig4a_part6 fig4b_part6


fig4_alphapower:
	-mkdir data/fig4
	python exp/fig4p.py data/fig4/4p pars/fig4/fig4p.yaml

# -

fig4a_part1:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part1 \
		pars/fig4/mathewson_constant_osc.yaml \
		-t 1

fig4a_part2:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part2 \
		pars/fig4/mathewson_constant_osc.yaml \
		-t 1.5

fig4a_part3:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part3 \
		pars/fig4/mathewson_constant_osc.yaml \
		-t 2.0

fig4b_part1:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part1 \
		pars/fig4/mathewson_lockedburst_osc.yaml \
		-t 1

fig4b_part2:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part2 \
		pars/fig4/mathewson_lockedburst_osc.yaml \
		-t 1.5

fig4b_part3:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part3 \
		pars/fig4/mathewson_lockedburst_osc.yaml \
		-t 2.0

# -
fig4_noosc: fig4a_part4 fig4a_part5 fig4a_part6 fig4b_part4 fig4b_part5 fig4b_part6


fig4a_part4:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part4 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 1

fig4a_part5:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part5 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 1.5

fig4a_part6:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/a_part6 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 2.0

fig4b_part4:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part4 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 1

fig4b_part5:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part5 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 1.5

fig4b_part6:
	-mkdir data/fig4
	nice -n 19 python exp/fig4.py data/fig4/b_part6 \
		pars/fig4/mathewson_constant_noosc.yaml \
		-t 2.0

