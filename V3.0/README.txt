########################################################################################################################
########################################################################################################################
###############################################  ORTHOLOGY_UTILS  ######################################################
########################################################################################################################
########################################################################################################################
The following programs rely on the location of the Orthology_utils folder to be known by the python environment. This
means updating the PYTHONPATH. To do this :
# On windows :
Go to : Settings > Advance system Settings > Environment variables and add PYTHONPATH to the user variables with a path
for the Orthology_utils folder. Only one reference to PYTHONPATH must exist, if multiple packages, uses ";" to delimit
their location.
# On Linux :
Open the ~/.bashrc file, a hidden file in the home folder of the user. Add the line :
export PYTHONPATH=$PYTHONPATH:/path/to/Orthology_utils
# Alternatively, you can add the path to Orthology utils to the PYTHONPATH of your IDE

########################################################################################################################
########################################################################################################################
############################################  PROGRAM DESCRIPTIONS  ####################################################
########################################################################################################################
########################################################################################################################

These programs use command line parsers to specify execution options. Each option can be specified using either the long
 (starts with two dashes) or short flag for the options. Flags can be necessary (no default), or optional, in which case
 a default value is used. In the case of input optional file name, if no name is provided, a file with the default name
  must exist in the parent directory. When specified, multiple option can be given for a single flag. All specie names
  must respect the naming convention in gProfiler (i.e. the id on https://biit.cs.ut.ee/gprofiler/page/organism-list )

#############################                       1_Get_orthologs                        #############################
######## Description
This programs takes a list of genes for a given species, and outputs the orthologs for other species, as well as
retrieves the sequences associated with each pair of orthologs.

######## Arguments
#### Required
## Input
-f   --gene_file           Name of test file containing the name of the genes in the "original" species in the first
                           column, other columns are ignored.
-os  --original_specie     Name of the specie corresponding to the gene names provided in the input file.
-as  --additional_species  Additional species. Can take multiple arguments.

#### Optional
## Output
-of  --orthology_file      Name of the orthology database json output file. Default is Orthology.json
-sf  --sequences_file      Name of the sequences database json output file. Default is Sequences.json
-nf  --names_file          Name of the names database json output file. Default is Names.json

############################################  2_Get_IDRs  ####################################################
######## Description
The programs looks at each individual sequences for orthologs, and predicts the disordered and folded domains location
from the primary protein sequence.

######## Arguments
#### Optional
##Input
-sf  --sequences_file      Name of the sequences database json output file. Default is Sequences.json
##Output
-pf  --properties_file     Name of the sequence properties file. Default is Sequence_properties.json.

############################################  3_Add_sequence_properties  ###############################################
######## Description
Adds sequence properties to the property file. Sequence properties can be : (1) sequence features (fraction of charge
residue (FCR), fraction of positive and negative residues, Net Charge Per Residue from primary sequence (Does not
account for pH)) (2) Sequence ensemble properties prediction from Sparrow (only strictly valid for IDRs), which include
end to end distance, radius of gyration, asphericity, ect... or (3) pH specific charge properties prediction based on
MEDOC (Isoelectric point and Net Charge, although mostly valid for disordered regions, is always done for everything).

######## Arguments
#### Optional
-of  --orthology_file      Name of the orthology database json output file. Default is Orthology.json
-sf  --sequences_file      Name of the sequences database json output file. Default is Sequences.json
-pf  --properties_file     Name of the sequence properties file. Default is Sequence_properties.json
-dsf --do_seq_feat         Whether to analyze sequence features. 0 if off, 1 is on. Default is 1 (on)
-dse --do_seq_ens          Whether to predict sequence ensemble properties. 0 if off, 1 is on. Default is 1 (on)
-dsp --do_seq_pH           Whether to perform pH calculations. 0 if off, 1 is on. Default is 1 (on)
-fp  --force_pred          Whether to force ensemble calculations on every domain types. By default, the ensemble
                           prediction are limited to all non-folded domains. 0 if off, 1 is on. Default is 0 (off)
-pHv --pH_val              Values at which to calculate the effective charge. Can have multiple arguments.
###############################################  4_Compare_sequences  ##################################################
######## Description
Compares sequences of orthologs between a reference species and a number of other species. Comparisons includes
comparing sequence features, sequence ensemble prediction and charge sequence features, as well as calculating the
homology between sequences.

######## Arguments
#### Required
-rs  --reference_specie    Name of the reference specie to which sequences are compared
-as  --additional_species  Name of additional species which are compared to the reference specie. Can have multiple
                           arguments.

#### Optional
-of  --orthology_file      Name of the orthology database json output file. Default is Orthology.json
-sf  --sequences_file      Name of the sequences database json output file. Default is Sequences.json
-pf  --properties_file     Name of the sequence properties file. Default is Sequence_properties.json
-hf  --homologies_file     Name of the homology database file. Default is Gene_homologies.json
-mlf --min_len_fraction    Minimum length fraction. Best to keep 0, build the database for everything, and discard later
                           , however, if the algorithm runs too slowly, this is a simple way of speeding it up.
-dh  --do_homo             Whether to calculate the homology score. 1 is on, 0 is off. Default is 1 (on).
-at  --align_type          Alignment type for the homology. Can be NW (Needleman–Wunsch, i.e. global) or SW
                           (Smith–Waterman, i.e. local). Default is SW.

################################################  5_Plots_comparison ###################################################
######## Description
Plot the sequence comparison between two species with regard to a third reference specie.
######## Arguments
#### Required
-rs  --reference_specie    Name of the reference specie to which sequences are compared
-ts  --top_specie          Top specie for the comparison (top_specie/norm_specie)
-ns  --norm_specie         Norm specie for the comparison (top_specie/norm_specie)
#### Optional
-hf  --homologies_file     Name of the homology database file. Default is Gene_homologies.json
-nf  --names_file          Name of the names database json output file. Default is Names.json
-tif --top_iso_fraction    Top fraction of isoforms that are kept using overall homology as a metric. 0 is only most
                           homologous, 1 is all.
-tof --top_orth_fraction   Top fraction of orthologs that are kept using overall homology as a metric. 0 is only most
                           homologous, 1 is all.
-mlr --min_len_ratio       Minimum length ratio between the species and the reference species for orthologs to be
                           plotted
-bw  --bin_width           Width of the bins for histograms. Default is 0.2
-flr --factor_len_ratio    Whether to multiply the homology by the length ratio. Can be 1 (on) or 0 (off). Default is 0

################################################  6_Plot_properties ####################################################
######## Description
This program reads the Orthology and Sequence properties files and plot all relevant properties individually for each
species.
######## Arguments
#### Required
-of  --orthology_file      Name of the orthology database json output file. Default is Orthology.json
-pf  --properties_file     Name of the sequence properties file. Default is Sequence_properties.json

################################################  7_Check_existence ####################################################
######## Description
This program takes the orhtology database and checks which species have an orthologs
######## Arguments
#### Required
-s  --species      Names of the species to be checked.

########################################################################################################################
########################################################################################################################
###################################################  DATABASES  ########################################################
########################################################################################################################
########################################################################################################################
There are two distinct databases that the algorithm works with :
    -The orthology database : contains data about the orthologs and name of each genes
    -The sequence database : contains the sequence for each isoform ID
    -The homology database : contains all data about sequence comparisons between homologs
    -The names database : contains the gene names for each gene ID

#############ORTHOLOGY DATABASE#############
Structure :
Original specie gene name > Specie > Gene ID > Protein isoform ID

#############SEQUENCE DATABASE#############
Structure :
Protein isoform ID > Amino acid sequence

#############HOMOLOGY DATABASE#############
reference specie > Compared specie > reference ortholog ID > Compared ortholog ID > reference protein isoform ID >
compared protein ioform ID > domain annotation > compared region pair > compared value type > value

#############NAMES DATABASE#############
Gene ID > Gene name

########################################################################################################################
########################################################################################################################
###################################################  DEPENDENCIES  #####################################################
########################################################################################################################
########################################################################################################################
Necessary python3 modules for this are :
Import name     Full name       Link
argparse        argparse        https://pypi.org/project/argparse/
requests        requests        https://pypi.org/project/requests/
numpy           numpy           https://pypi.org/project/numpy/
metapredict     metapredict     https://pypi.org/project/metapredict/
matplotlib      matplotlib      https://pypi.org/project/matplotlib/
skbio           scikit-bio      https://pypi.org/project/scikit-bio/
sparrow         sparrow         https://sparrow-online.com/

To install all on linux type :

pip install -r requirements.txt

########################################################################################################################
########################################################################################################################
###############################################  EXECUTION EXAMPLE  ####################################################
########################################################################################################################
########################################################################################################################
python3 1_Get_orthologs_V1.3.py -os drerio -as mmusculus hsapiens -f batch_39.txt
python3 2_Get_IDRs_V2.0.py
python3 3_Add_sequence_properties_V1.1.py
python3 4_Compare_sequences_V1.3.py -rs hsapiens -as mmusculus drerio
python3 5_Plot_comparison_V1.0.py -rs hsapiens -ts drerio -ns mmusculus
python3 6_Plot_properties_V2.4.py

An example of batches execution for a slurm submission system is given in Batches_Setup/Setup.sh
You can use the Recombine_batches_json_V1.0.py -n with the nax number of batches to recombine the batches into a single
file per database type
