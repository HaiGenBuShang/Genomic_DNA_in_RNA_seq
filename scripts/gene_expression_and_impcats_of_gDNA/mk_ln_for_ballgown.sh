stringtie_ballgown_dir=${1}
ballgown_for_R_dir=${2}
old_and_new_dir_table=${3}

while read i j
do
ln -s ${stringtie_ballgown_dir}/${i} ${ballgown_for_R_dir}/${j}
done < <(cat ${old_and_new_dir_table})
