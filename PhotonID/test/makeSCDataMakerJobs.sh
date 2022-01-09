#./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
EXECUTABLE=mvaDataMaker.exe
CFG_TEMPLATE=configs/QcdSample.template
ANALYSIS_OPT=1

declare -a SourceFiles=(\
"srcFiles/QCD_Pt_15to30_TuneCP5_13TeV_pythia8.files" \
"srcFiles/QCD_Pt_30to50_TuneCP5_13TeV_pythia8.files" \
)

declare -a tagArr=(\
"mc_QCD_Pt_15to30" \
"mc_QCD_Pt_30to50" \
)

declare -a AnalysisOption=(\
1 \
1 \
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
        200 \
        2 \
        $TAG
#    set +x
done
