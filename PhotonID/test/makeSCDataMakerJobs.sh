#./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
NJOBS=${1-5}
FILES_PER_JOB=${2-1}
MAXEVENTS=${3--1}
echo NJOBS : $NJOBS
echo FILES_PER_JOB : $FILES_PER_JOB
echo MAXEVENTS : $MAXEVENTS
echo ""
EXECUTABLE=mvaDataMaker.exe
#CFG_TEMPLATE=configs/QcdSample_template.cfg

declare -a SourceFiles=(\
"srcFiles/bs2mmg.files" \
"srcFiles/flatpi0.fls" \
"srcFiles/qcd15To7000.fls" \
"srcFiles/qcd20To30EMEnriched.fls" \
"srcFiles/qcd30To50EMEnriched.fls" \
"srcFiles/qcd30To50.fls" \
"srcFiles/qcd50To80.fls" \
"srcFiles/data_2018D.fls" \
)

declare -a tagArr=(\
"mc_sig_bs2mmg" \
"mc_bkg_flatpi0" \
"mc_bkg_qcd15To7000" \
"mc_bkg_qcd20To30EMEnriched" \
"mc_bkg_qcd30To50EMEnriched" \
"mc_bkg_qcd30To50" \
"mc_bkg_qcd50To80" \
"data_2018D" \
)

declare -a AnalysisOption=(\
2 \
1 \
1 \
1 \
1 \
1 \
1 \
3 \
)

declare -a CfgTemplate=(\
"configs/BsToMMGSample.tpl" \
"configs/QcdSample.tpl" \
"configs/QcdSample.tpl" \
"configs/QcdSample.tpl" \
"configs/QcdSample.tpl" \
"configs/QcdSample.tpl" \
"configs/QcdSample.tpl" \
"configs/data.tpl" \
)

# ./makeCondorJobForAnalysis.py mvaDataMaker.exe configs/QcdSample.tpl 1 /home/athachay/t3store3/bs2mumug/photonID/analysis/CMSSW_10_6_29/src/BsMMGAnalysis/PhotonID/test/results/MC/mc_bkg_qcd30To50EMEnriched 5 1 -1 mc_bkg_qcd30To50EMEnriched
for i in "${!tagArr[@]}"; do 
    echo $i : ${jobArr[$i]}
#    set -x
    src=${SourceFiles[$i]}
    TAG=${tagArr[$i]}
    ANALYSIS_OPT=${AnalysisOption[$i]}
    CFG_TEMPLATE=${CfgTemplate[$i]}
#    set +x
    #set -x
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
    #set +x
done
