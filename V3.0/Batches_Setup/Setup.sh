source ~/.bashrc
batch_size=25
file_name=Gene_names.txt
SCRIPT_LOC='../'

mkdir Batches
mkdir Batches/Input_names
mkdir Batches/Names
mkdir Batches/Sequences
mkdir Batches/Ortho
mkdir Batches/Seq_prop
mkdir Batches/Homology
mkdir Batches/Standard_output
awk -v size="${batch_size}" '{print $1 > ("./Batches/Input_names/batch_" int((NR-1)/size) ".txt")}' ${file_name}
N_batches=$(( ($(wc -l < "$file_name") + batch_size - 1) / batch_size ))

for((i=0;i<${N_batches};i++));
do
seq_file=./Batches/Sequences/N_${i}.json
ortho_file=./Batches/Ortho/N_${i}.json
names_file=./Batches/Names/N_${i}.json
homo_file=./Batches/Homology/N_${i}.json
seq_prop_file=./Batches/Seq_prop/N_${i}.json

echo "#!/bin/bash" > run.batch
echo "#SBATCH -o Batches/Standard_output/std"$i".out" >> run.batch
echo "#SBATCH -e Batches/Standard_output/std"$i".err" >> run.batch
echo "#SBATCH --ntasks=1" >> run.batch
echo "#SBATCH --cpus-per-task=1"  >> run.batch

echo python3 ${SCRIPT_LOC}+/1_Get_orthologs_V1.3.py -f ./Batches/Input_names/batch_${i}.txt -of ${ortho_file} -sf ${seq_file} -nf ${names_file} -os drerio -as hsapiens mmusculus>> run.batch
echo python3 ${SCRIPT_LOC}+/2_Get_IDRs_V2.0.py -of ${ortho_file} -sf ${seq_file} -pf ${seq_prop_file} >> run.batch
echo python3 ${SCRIPT_LOC}+/3_Add_sequence_properties_V1.1.py -pf ${seq_prop_file} -sf ${seq_file} -of ${ortho_file}  >> run.batch
echo python3 ${SCRIPT_LOC}+/4_Compare_sequences_V1.3.py -of ${ortho_file} -hf ${homo_file} -sf ${seq_file} -pf ${seq_prop_file} -rs hsapiens -as drerio mmusculus >> run.batch
#echo python3 ${SCRIPT_LOC}+/5_Plot_comparison_V1.0.py -nf ${names_file} -hf ${homo_file} -rs hsapiens -as drerio mmusculus >> run.batch
#echo python3 ${SCRIPT_LOC}+/6_Plot_properties_V2.4.py -f ./Batches/Homology_scores/N_${i}_All_homology.json -rs hsapiens  -as drerio mmusculus >> run.batch

sbatch -J N_${i} run.batch
sleep 1
done
