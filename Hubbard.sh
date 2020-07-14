#!/bin/dash

chmod +x py_scripts/Hubbard.py

#mode : {fromNI, fromIC, refine}

FLAGS="--cmix=0.8 --Niter=50 --Msteps=5e5 --NC=52 --mode=refine --nom=80 --adiabatic=True "

TWO=2.0

#Set up spooler

#Set up logging

mkdir -p log

#Metals

prefix="Hubbard/Metal"

for T in 0.03;
do
	rm initialcondition/*
	rm run/*
	#cp results/$prefix/T.0.0300/U.2.2000/Gf.out initialcondition/Gf.out
	prefix="Hubbard.test/Metal"
	for U in 2.5;
	do
		mu=$(echo "scale=4 ; $U / $TWO" | bc)
		current_time=$(date "+%Y.%m.%d-%H.%M.%S") #timestamp for logfile
		SPATH=$(printf "$prefix/T.%.04f/U.%.04f/" $T $U)
		#py_scripts/Hubbard.py --mu=$mu --T=$T --U=$U --spath=$SPATH $FLAGS
		py_scripts/Hubbard.py --mu=$mu --T=$T --U=$U --spath=$SPATH $FLAGS | tee log/log.$current_time
	done
done

prefix="Hubbard/Insul"

for T in 0.03;
do
	rm initialcondition/*
	rm run/*
	#cp results/$prefix/T.0.0300/U.4.0000/Gf.out initialcondition/Gf.out
	prefix="Hubbard.test/Insul"
	for U in 2.5;
	do
		mu=$(echo "scale=4 ; $U / $TWO" | bc)
		current_time=$(date "+%Y.%m.%d-%H.%M.%S") #timestamp for logfile
		SPATH=$(printf "$prefix/T.%.04f/U.%.04f/" $T $U)
		py_scripts/Hubbard.py --mu=$mu --T=$T --U=$U --spath=$SPATH $FLAGS | tee log/log.$current_time
	done
done



