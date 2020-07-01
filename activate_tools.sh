tool=$1

source /users/brc/whh839/.bashrc
#module load Anaconda3/5.3.0

if [ "$#" -ne 1 ] 
then
    echo "Illegal number of arguments
    use: activate_tools.sh tool_name 
    use: activate_tools.sh list to show a list of possible options"    
fi

if [ "$tool" == "CNV" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/CNV_analysis
elif [ "$tool" == "scikit" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/VCF_scikit
elif [ "$tool" == "delly" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/delly
elif [ "$tool" == "gatk" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/gatk
elif [ "$tool" == "mosdepth" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/moshdepth
elif [ "$tool" == "vcf2db" ]
then
    conda activate /well/gel/HICF2/software/anaconda3/envs/vcf2db
elif [ "$tool" == "CNV" ]
then
    conda activate /well/gel/HICF2/software/gemini/data/anaconda
elif [ "$tool" == "list" ]
then
    echo "List of possibile selections
    1. CNV: CNV related tools (SVtools, SVtyper, lumpy) 
    2. scikit: python module scikit-allel for VCF manipulation
    3. delly: delly CNV / SV caller
    4. gatk: GATK 4 with deep learning variant filter module
    5. mosdepth: mosdepth tools for coverage analysis
    6. vcf2db: vcf2db tool for loading VCF into GEMINI compatible db
    7. gemini: gemini variant dataset manager"
else
    echo "Unrecognized tool name
    use: activate_tools.sh list to show a list of possible options"
fi
