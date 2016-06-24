#!/bin/bash
if [ "$1" = "slo" ] || [ "$1" = "ido" ] || [ "$1" = "rlf" ] || [ "$1" = "lfo" ] || [ "$1" = "sdo" ]; then
	for file in data/*.mtx
	do
	 # do something on "$file"
		#file name is now like data/Matrix.mtx
		#the following lines are used to remove .mtx and /data	  
		#split file name using '.' delimeter 
		IFS='.'
		read -a direc <<< "${file}"
		#split file name using '.' delimeter 
		IFS='/'
		read -a name <<< "${direc}"

		#make necessary directories
		mkdir -p Profiling_Result
		mkdir -p Profiling_Result/slo
		mkdir -p Profiling_Result/ido
		mkdir -p Profiling_Result/rlf
		mkdir -p Profiling_Result/lfo
		mkdir -p Profiling_Result/sdo

		#Put gcolor result and profilings in appropriate folders
		../examples/gcolor -i "$file" -m "$1" > Profiling_Result/"$1"/${name[1]}_Result.out
		cd ../examples
		gprof gcolor > ../dsjm_testing/Profiling_Result/"$1"/${name[1]}_Prof.out
		cd ../dsjm_testing
	done
else
	echo "Wrong method! Please pass slo/ido/rlf/lfo/sdo as argument"
fi
	#echo $1
	
	
 

