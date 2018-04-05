SHELL=/bin/bash -O expand_aliases
# type
# export ALPHA_DATADIR=/Users/type/Data/Smith/data
# export ALPHA_CACHEDIR=/Users/type/Data/Smith/data/tmp
# FNB
export ALPHA_DATADIR=/home/ejp/data/Smith/data
export ALPHA_CACHEDIR=/home/ejp/data/Smith/data/tmp

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


# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# SMITH DATA ANALYSIS
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------
clean_smith_cache:
	-rm -rf $(ALPHA_CACHEDIR)

# -----------------------------------------------------------------------------------------------
smith_analysis1:
	# Create cache
	-mkdir $(ALPHA_CACHEDIR)
	# Clean old results
	-rm $(ALPHA_DATADIR)/analysis1_*.csv
	# Run
	parallel -j 12 -v \
		--joblog '$(ALPHA_DATADIR)/analysis1.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis1_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.02' ::: {0..20}
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)

# -----------------------------------------------------------------------------------------------
# 4-2-2018
# Made some cache opt since a1.
# Here's another 'quick' test.
# 88ecc44df7b2bb061796b33c0b7313a31323c716
# RUnning all sets at the time let the cache grow to > 300 Gb
# Going to need to run in smaller sets, clearing the cache manually along the way.
smith_analysis2:
	# Create cache
	-mkdir $(ALPHA_CACHEDIR)
	# Clean old results
	-rm $(ALPHA_DATADIR)/analysis2_*.csv
	# Run
	parallel -j 12 -v \
		--joblog '$(ALPHA_DATADIR)/analysis2.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis2_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.05' ::: {0..20}
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)

# -----------------------------------------------------------------------------------------------
# 4-3-2018
# Try it in smaller chunks..., 
# out of order to sample both monkeys in the in the first couple runs. 
# Sampling 0.2 of the total data each run will take 2-3 days to finish...
# 53cba2facfcb2e4fd4939067e4f25ed876adc189
smith_analysis3:
	# ------------------------------
	# Create cache
	-mkdir $(ALPHA_CACHEDIR)
	-rm -rf $(ALPHA_CACHEDIR)
	# Clean old results
	-rm $(ALPHA_DATADIR)/analysis3_*.csv
	# ------------------------------
	# Run 1 - 5
	parallel -j 5 -v \
		--joblog '$(ALPHA_DATADIR)/analysis3.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis3_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.2' ::: 0 1 2 4 5
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)
	# ------------------------------
	# Run 10 - 15
	parallel -j 5 -v \
		--joblog '$(ALPHA_DATADIR)/analysis3.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis3_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.2' ::: 11 12 13 14 15
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)
	# ------------------------------
	# Run 5 - 10
	parallel -j 5 -v \
		--joblog '$(ALPHA_DATADIR)/analysis3.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis3_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.2' ::: 6 7 8 9 10
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)
	# ------------------------------
	# Run 15 - 21
	parallel -j 5 -v \
		--joblog '$(ALPHA_DATADIR)/analysis3.log' \
		--nice 19 --delay 2 --colsep ',' \
	'python exp/smith_burst_analysis.py analysis3_set{1} $(ALPHA_DATADIR) --verbose --n={1} --percent_segment=0.2' ::: 16 17 18 19 20
	# Cleanup cache
	-rm -rf $(ALPHA_CACHEDIR)
	