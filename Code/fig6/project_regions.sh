#! /bin/bash

# This script takes a bed file of point coordinates and projects them from a reference to a query.
# Make sure path_pwaln_pkl contains the desired species set.

[[ $# < 6 ]] && { echo "Usage: ./project_regions.sh coords.bed reference_species query_species half_life_distance nthreads outdir [optional: path_pwaln_pkl]"; exit 0; }

bed=$1
ref=$2
qry=$3
half_life_distance=$4
nthreads=$5
outdir=$6

# fetch pairwise alignment pickle
if [[ $# == 7 ]]; then
	path_pwaln_pkl=$7
else
	path_pwaln_pkl=/project/ipp-data/projection/${ref}_${qry}/${ref}.${qry}.pwaln.pkl
	[[ ! -f $path_pwaln_pkl ]] && path_pwaln_pkl=/project/ipp-data/projection/${qry}_${ref}/${qry}.${ref}.pwaln.pkl
	[[ ! -f $path_pwaln_pkl ]] && echo -e "I did not find a collection of pairwise alignments for this species pair at /project/ipp-data/projection.\nCreate one with a given selection of bridging species using preprocess_wrapper.sh" && exit 1
fi

# check if bed-file contains at least 4 columns and unique names in column 4
[[ $(awk '{print NF}' $bed | sort -nu | head -n 1) -lt 4 ]] && echo "Error: bed-file must contain at least 4 columns with the 4th being a unique ID / name." && exit 1
[[ $(cut -d$'\t' -f4 "$bed" | sort | uniq | wc -w) != $(< "$bed" wc -l) ]] && echo "Error: Names in column 4 of bed-file must be unique." && exit 1

tmp_dir=$(sed -e 's/bed/tmp/' -e 's/\.bed//' <<< "$bed")
mkdir -p $tmp_dir
l=$(< $bed wc -l)
ERT=$(printf "%.0f" "$(echo "1*$l/$nthreads" | bc)") # based on a estimated average runtime of 1 min per job
echo "Estimated runtime: $((ERT/60))h $(bc <<< $ERT%60)m"
sem_id="project_coordinates_$(hostname)_${RANDOM}"
starttime=$(date -u '+%s')

# loop through bed-file
i=0
N=$(wc -l $bed)
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
	bed_row=($LINE)
	id=${bed_row[3]}
	# ((i++)) # implement progress bar instead of printing every single coordinate
	# f=$(bc <<< 'scale=2; $i/$N')
	### continue here printf '=%.0s' {1..100}
	# echo -ne '($(i/N)%)\r' 
	[[ -f ${tmp_dir}/${id}.proj.tmp ]] && echo "${tmp_dir}/${id}.proj.tmp exists. Skip." && continue
	coord=${bed_row[0]}:$(($((${bed_row[1]}+${bed_row[2]}))/2)) # center of region
	echo $id $coord
	sem --id $sem_id -j${nthreads} project_dijkstra.py $ref $qry $coord $id $half_life_distance $path_pwaln_pkl $tmp_dir # sem is an alias for parallel --semaphore. A counting semaphore will allow a given number of jobs to be started in the background.
done < $bed

sem --id $sem_id --wait # wait until all sem jobs are completed before moving on
endtime=$(date -u '+%s')
difftime=$(date -u --date @$((endtime-starttime)) '+%-Hh %-Mm %-Ss')
echo "Effective runtime: ${difftime}"

# concatenate tmp output files to one file, delete tmp files
outfile=${outdir}/$(basename $bed .bed).proj
mkdir -p $(dirname $outfile)
echo "Concatenating temporary files and writing to ${outfile}"
ids=($(cut -f4 $bed))
head -3 "${tmp_dir}/${ids[0]}.proj.tmp" > $outfile
# echo -n "" > ${bed}_unmapped
for id in ${ids[@]:1}; do
	if [ -f "${tmp_dir}/${id}.proj.tmp" ]; then
		eval tail -n 1 -q "${tmp_dir}/${id}.proj.tmp" >> $outfile
	else
		echo $id >> ${outfile}_unmapped
	fi
done
rm -r $tmp_dir

echo "Done"
