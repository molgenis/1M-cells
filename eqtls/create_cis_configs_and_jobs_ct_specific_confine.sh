#!/usr/bin/env bash
FEATURES_LOC=$1
FEATURES_LOC_2=$2
CONFIGS_LOC=$3
RESULTS_LOC=$4
JOBS_LOC=$5
CONFINEMENTS_PREPEND=$6
CONFINEMENTS_APPEND=$7

cd ${FEATURES_LOC}
for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done

# loop through files
files=${FEATURES_LOC}"/*"
for f in ${files}
do
  # only create files for features
  if [[ ${f} =~ .*\.(tsv) ]]
  then
    # ${f} is the features file full path, we also need just the name
    # echo "${f}"
    fbasename=$(basename "${f}" .tsv)
    expression_files_loc_2=${FEATURES_LOC_2}"/"${fbasename}".tsv"
    # echo "${fbasename}"
    # make basename posix compliant (i.e. no spaces etc.)
    posix_basename=${fbasename//[^a-zA-Z0-9]/_}
    #posix_basename="echo ${fbasename} | tr ' ' '_'"
    #mv $f "${FEATURES_LOC}/${posix_basename}.tsv"
    # echo ${posix_basename}
    # create location for config file
    config_location=${CONFIGS_LOC}"/"${posix_basename}"_emp_config.xml"
    # create the folder if non_existant
    mkdir -p ${CONFIGS_LOC}
    echo ${config_location}
    # create location for the result files
    result_location=${RESULTS_LOC}"/"${posix_basename}"/"
    # create the folder if non_existant
    mkdir -p ${result_location}
    echo ${result_location}
    # create location for job file
    job_location=${JOBS_LOC}"/"${posix_basename}"_emp_job_SBATCH.sh"
    # create the folder if non_existant
    mkdir -p ${JOBS_LOC}
    echo ${job_location}
    # the confinements don't have '_expression' in them, so we need to strip that
    expressionless_ct=$(echo ${fbasename} | sed 's/_expression//g'

  # create the config file
echo "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>
<settings>
  <defaults>
    <qc>
        <snpqccallratethreshold>0.95</snpqccallratethreshold>
        <snpqchwethreshold>0.0001</snpqchwethreshold>
        <snpqcmafthreshold>0.1</snpqcmafthreshold>
    </qc>
    <analysis>
        <analysistype>cis</analysistype>
        <cisanalysisprobedistance>100000</cisanalysisprobedistance>
        <correlationtype>nonparametric</correlationtype>
        <!-- <regressOutEQTLEffects>/groups/umcg-bios/scr01/umcg-aclaringbould/eqtlpipeline_lld/cis-eQTLs/conditional/output_round18_RegressedOut/toRegressOut.txt</regressOutEQTLEffects> -->
        <threads>10</threads>
        <createdotplot>false</createdotplot>
        <createqqplot>false</createqqplot>
    </analysis>
    <multipletesting>
        <type>fdr</type>
        <threshold>0.05</threshold>
         <fdrtype>probe-level</fdrtype>
        <permutations>10</permutations>
    </multipletesting>
    <output>
        <outputdirectory>${result_location}</outputdirectory>
        <outputplotthreshold>0</outputplotthreshold>
        <outputplotdirectory>${result_location}</outputplotdirectory>
        <maxnreqtlresults>50000000</maxnreqtlresults>
        <generatesnpsummarystatistics>false</generatesnpsummarystatistics>
        <generateeqtlpvaluetable>false</generateeqtlpvaluetable>
        <binaryoutput>false</binaryoutput>
        <textoutput>true</textoutput>
    </output>
    <confine>
        <snp/>
	<snpProbe>${CONFINEMENTS_PREPEND}${expressionless_ct}${CONFINEMENTS_APPEND}</snpProbe>
        <confineSNPsToSNPsPresentInAllDatasets>true</confineSNPsToSNPsPresentInAllDatasets>
        <confineSNPsSelectSNPInStrongestLD>false</confineSNPsSelectSNPInStrongestLD>
        <confineProbesThatMapToKnownChromosome>true</confineProbesThatMapToKnownChromosome>
    </confine>
  </defaults>
  <datasets>
    <dataset>
        <name>${posix_basename}_v2</name>
		<location>/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/genotypes/LL_tritypes/</location>
		<expressiondata>${f}</expressiondata>
		<probeannotation>/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/eQTLgen-rebuttal/trans_eqtls/singleCell-annotation-stripped.tsv</probeannotation>
		<quantilenormalize>false</quantilenormalize>
		<logtranform>false</logtranform>
    </dataset>
    <dataset>
        <name>${posix_basename}_v3</name>
        <location>/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/genotypes/LL_tritypes/</location>
        <expressiondata>${expression_files_loc_2}</expressiondata>
        <probeannotation>/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/eQTLgen-rebuttal/trans_eqtls/singleCell-annotation-stripped.tsv</probeannotation>
        <quantilenormalize>false</quantilenormalize>
        <logtranform>false</logtranform>
    </dataset>
  </datasets>
</settings>" > ${config_location}

# create the job file
echo "#!/usr/bin/env bash
#SBATCH --job-name=map_${posix_basename}
#SBATCH --output=${JOBS_LOC}/map_${posix_basename}.out
#SBATCH --error=${JOBS_LOC}/map_${posix_basename}.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
set -e
ml Java
java -jar -Xmx40g -Xms20g -XX:StringTableSize=10000019 -XX:MaxPermSize=512m \
/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/eQTLgen-rebuttal/trans_eqtls/eqtl-mapping-pipeline-1.4.0-SNAPSHOT/eqtl-mapping-pipeline.jar \
--mode metaqtl \
--settings ${config_location}" > ${job_location}
  fi
done
