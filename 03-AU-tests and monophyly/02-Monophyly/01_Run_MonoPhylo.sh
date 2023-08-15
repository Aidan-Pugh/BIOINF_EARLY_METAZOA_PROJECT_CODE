#!/bin/bash

# Iterate over all directories starting with "tree_"

for directory in tree_*; do
	if [ -d "$directory" ]; then
		cd "$directory"
		log_file_times="log_file_times.txt"
		touch $log_file_times

		echo "Started at: " >> "$log_file_times"
		date +%H:%M:%S >> "$log_file_times"
		echo "----------------------- Found $directory: Running Monophylo ----------------------------------"
		
		file=$(find . -maxdepth 1 -type f -name "renamed_*")


		python3 ../MonoPhylo.py -t /home/your/pathway/directory/Monophylo/"$directory"/$file -o /home/your/pathway/directory/Monophylo/"$directory"/Results/ -m /home/your/pathway/directory/Monophylo/"$directory"/taxonomy_file.txt

		echo "Monophylo complete in $directory"
		echo "Finished at: " >> "$log_file_times"
		date +%H:%M:%S >> "$log_file_times"

		cd ..
	fi
done
