#!/bin/dash

chmod +x py_scripts/Hubbard.py
mkdir -p log

#mode : {fromNI, fromIC, refine}

FLAGS="--cmix=0.8 --Niter=20 --Msteps=4e6 --NC=50 --mode=refine --nom=500 --adiabatic=True --U=4.0"

prefix="DopedHubbard.Lostat/Insul"

for T in 0.005 0.004 0.003 0.002;
do
	rm initialcondition/*
	rm run/*
	for mu in 2.0 2.2 2.4 2.6 2.8 2.9 2.91 2.92 2.93 2.94 2.95 2.96 2.97 2.98 2.99 3.0;
	do
		current_time=$(date "+%Y.%m.%d-%H.%M.%S") #timestamp for logfile
		SPATH=$(printf "$prefix/T.%.04f/mu.%.04f/" $T $mu)
		py_scripts/Hubbard.py --mu=$mu --T=$T --spath=$SPATH $FLAGS | tee log/log.$current_time
	done
done

FLAGS="--cmix=0.5 --Niter=20 --Msteps=4e6 --NC=50 --mode=fromIC --nom=500 --adiabatic=True --U=4.0"

prefix="DopedHubbard.Lostat/Metal"

for T in 0.005 0.004 0.003 0.002;
do
	rm initialcondition/*
	rm run/*
	for mu in 3.2 3.1 3.0 2.99 2.98 2.97 2.96 2.95 2.94 2.93 2.92 2.91 2.9 2.8;
	do
		current_time=$(date "+%Y.%m.%d-%H.%M.%S") #timestamp for logfile
		SPATH=$(printf "$prefix/T.%.04f/mu.%.04f/" $T $mu)
		py_scripts/Hubbard.py --mu=$mu --T=$T --spath=$SPATH $FLAGS | tee log/log.$current_time
	done
done

FLAGS="--cmix=0.8 --Niter=20 --Msteps=4e6 --NC=50 --mode=refine --nom=500 --adiabatic=True --U=4.0"

prefix="DopedHubbard.Lostat/Metal"

for T in 0.005 0.004 0.003 0.002;
do
	rm initialcondition/*
	rm run/*
	for mu in 3.2 3.1 3.0 2.99 2.98 2.97 2.96 2.95 2.94 2.93 2.92 2.91 2.9 2.8;
	do
		current_time=$(date "+%Y.%m.%d-%H.%M.%S") #timestamp for logfile
		SPATH=$(printf "$prefix/T.%.04f/mu.%.04f/" $T $mu)
		py_scripts/Hubbard.py --mu=$mu --T=$T --spath=$SPATH $FLAGS | tee log/log.$current_time
	done
done
