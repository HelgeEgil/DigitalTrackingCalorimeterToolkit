#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <file_number>"
echo File Number: $1

if [ $# -ne 1 ]; then
	echo Invalid number of arguments: $#
	exit
fi

nice -n 10 Gate -a "'[nfile,$1]" Main.mac > terminal_output.txt &
echo "Running: $1"
