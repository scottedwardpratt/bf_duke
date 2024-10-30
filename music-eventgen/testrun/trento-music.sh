#!/bin/bash
case $# in
0|1|2|3|4|5)
	echo "Usage: trento-music.sh  run_number qualifier b eta_kink trento_norm trento_p";
	exit 1 ;;
*)

################
## GRIDPARAMS ##
################
GRIDMAX=25
GRIDSTEP="0.200"
run_number=$1
qualifier=$2
b=$3
eta_kink=$4
trento_norm=$5
trento_p=$6

eta_lowT_slope=0.0
eta_highT_slope=0.0
bulk_norm=0.0
bulk_width=0.1
bulk_Tpeak=0.17

NEVENTS=1000

cat > ./trento_data/input-trento.dat << END

# specify the projectile option twice
projectile = Pb
projectile = Pb
number-events = $NEVENTS

# don't print event properties to stdout, save to HDF5
quiet = true
output = PbPb.hdf

reduced-thickness = ${trento_p}
fluctuation = 1.5
nucleon-width = 0.6
cross-section = 6.4
normalization = ${trento_norm}

# leave commented out for min-bias
b-min = ${b}
b-max = ${b}

grid-max = $GRIDMAX
grid-step = $GRIDSTEP
END

#echo Trento started running
rm -r ./trento_data/to_average/*
../models/trento/trento -c ./trento_data/input-trento.dat -o ./trento_data/to_average

echo Trento finished running

GRIDPOINTS=250

./trento-average.py ./trento_data/to_average/ $NEVENTS $GRIDPOINTS $GRIDSTEP
wait
cp ./trento_data/music-input-from-trento.dat ./


#cd music-hydro
cat > PCE_inputfile.dat << END


###################################
# parameters to play with
###################################
#initial distribution
Initial_Distribution_input_filename music-input-from-trento.dat
#chemical equilibrium params
Initial_light_quark_fugacity 1  # initial fugacities 
Initial_strange_quark_fugacity 1
Light_quark_equilibration_time 0.5 # equilibration timescales (fm/c)
Strange_quark_equilibration_time 0.5

# transport coefficients
Viscosity_Flag_Yes_1_No_0 1     # turn on viscosity in the evolution
Include_Shear_Visc_Yes_1_No_0 1 # include shear viscous effect
Shear_to_S_ratio 0.1            # value of \eta/s
T_dependent_Shear_to_S_ratio 3  # flag to use temperature dep. \eta/s(T)
eta_over_s_T_kink_in_GeV 0.17
eta_over_s_low_T_slope_in_GeV ${eta_lowT_slope}
eta_over_s_high_T_slope_in_GeV ${eta_highT_slope}
eta_over_s_at_kink ${eta_kink}

Include_Bulk_Visc_Yes_1_No_0 1  # include bulk viscous effect
T_dependent_Bulk_to_S_ratio 2
bulk_viscosity_normalisation ${bulk_norm}
bulk_viscosity_width_in_GeV ${bulk_width}
bulk_viscosity_peak_in_GeV ${bulk_Tpeak}

s_factor 60                   # normalization factor read in
END
cat PCE_fixed.txt >> PCE_inputfile.dat

echo 'run_number='${run_number}
echo 'qualifier='${qualifier}
echo 'bmin='${b}
echo 'bmax='${b}
echo 'bulk_norm='${bulk_norm}
echo 'bulk_width='${bulk_width}
echo 'bulk_Tpeak='${bulk_Tpeak}
echo 'eta_kink='${eta_kin}
echo 'deta/dT_low='${eta_lowT_slope}
echo 'deta/dT_high='${eta_highT_slope}
echo 'trento_norm='${trento_norm}
echo 'trento_p='${trento_p}

time ../models/music-hydro/MUSIChydro ./PCE_inputfile.dat
mkdir -p hydrodata/modelruns/run${run_number}/${qualifier};

mv evolution_xyeta.dat hydrodata/modelruns/run${run_number}/${qualifier}/hydrodata.txt

esac