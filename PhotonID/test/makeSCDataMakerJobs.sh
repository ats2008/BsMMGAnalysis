#./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
NJOBS=${1-5}
FILES_PER_JOB=${2-1}
MAXEVENTS=${3--1}
echo NJOBS : $NJOBS
echo FILES_PER_JOB : $FILES_PER_JOB
echo MAXEVENTS : $MAXEVENTS
echo ""
EXECUTABLE=mvaDataMaker.exe
CFG_TEMPLATE=configs/QcdSample_template.cfg

declare -a SourceFiles=(\
"srcFiles/bs2mmg.files" \
"srcFiles/qcd30to50Extended.files" \
#"srcFiles/QCD_Pt_15to30_TuneCP5_13TeV_pythia8.files" \
#"srcFiles/QCD_Pt_30to50_TuneCP5_13TeV_pythia8.files" \
)

declare -a tagArr=(\
"mc_sig_bs2mmg" \
"mc_QCD_Pt_30to50ext" \
#"mc_QCD_Pt_15to30" \
#"mc_QCD_Pt_30to50" \
)

declare -a AnalysisOption=(\
2 \
1 \
#1 \
#1 \
)

for i in "${!tagArr[@]}"; do 
    echo $i : ${jobArr[$i]}
    src=${SourceFiles[$i]}
    TAG=${tagArr[$i]}
    ANALYSIS_OPT=${AnalysisOption[$i]}
 #   set -x
     echo ./makeCondorJobForAnalysis.py \
        $EXECUTABLE \
        $src \
        $CFG_TEMPLATE \
        $ANALYSIS_OPT \
        $PWD/results/MC/$TAG \
        $NJOBS \
        $FILES_PER_JOB \
        $MAXEVENTS \
        $TAG
#    set +x
done
