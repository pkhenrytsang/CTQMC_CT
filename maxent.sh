#!/bin/bash

prefix="results/CT"

for T in 0.01;
do
	for t in 0.35;
	do
		SPATH=$(printf "$prefix/T.%.04f/Metal/t.%.04f/converged/" $T $t)
		CPATH=$(pwd)
		CPPATH=$(pwd)/maxent_params.dat
		MPPATH=$(pwd)/model.dat
		cd $SPATH
		mkdir maxent
		cd maxent
		
		echo "entering maxent directory"
		
		cp $CPPATH .
		cp $MPPATH .
		
		saverage.py --nexecute ../Sig.out.5[5-9]
		
		maxent_run.py sig.inpx
		
		#cp Sig.out ../Sig.ret.out

		cd $CPATH
		
		echo "going back to working directory"
	done
done
