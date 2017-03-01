

########################################
# Determine what to do
########################################

do_T2WS=false
do_bestfit=false
do_generateToys=false
do_covMat=false


if [ $# -eq 0 ]; then
    echo "No additional arguments supplied"
    return

elif [ "$1" == "all" ]; then
    do_T2WS=true
    do_bestfit=true
    do_generateToys=true
    do_covMat=true

elif [ "$1" == "onlyT2WS" ]; then
    do_T2WS=true
elif [ "$1" == "onlybestfit" ]; then
    do_bestfit=true
elif [ "$1" == "onlygenerateToys" ]; then
    do_generateToys=true
elif [ "$1" == "onlycovMat" ]; then
    do_covMat=true

elif [ "$1" == "skip" ]; then
    do_generateToys=true
    do_covMat=true

else
    echo "Pass an option (all or skip)"
    return

fi


function chapter {
    echo
    echo "----------------------------"
    echo "---------------------------- $1"
    echo "----------------------------"
    echo 
}

########################################
# Specify what to do
########################################

NTOYS=-1

NBINS=7
NCATS=3

do_pT=true
do_recoFit=false
do_regularization=true


# Select the appropriate datacard

# if [ "$do_regularization" = false ] ; then
#     if [ "$do_pT" = true ]; then
#         DC="Datacard_13TeV_differential_pT_moriond17.txt"
#     else
#         DC="Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod.txt"
#     fi
# else
#     if [ "$do_pT" = true ]; then
#         DC="Datacard_13TeV_differential_pT_moriond17_regularization.txt"
#     else
#         DC="Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_regularized.txt"
#     fi
# fi


if [ "$do_pT" = true ]; then

    if [ "$do_regularization" = false ] ; then
        DC="Datacard_13TeV_differential_pT_moriond17.txt"
    else
        DC="Datacard_13TeV_differential_pT_moriond17_regularization.txt"
    fi

    SQRT_TAU=0.035

else

    if [ "$do_regularization" = false ] ; then
        DC="Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod.txt"
    else
        DC="Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_regularized.txt"
    fi

    SQRT_TAU=0.78

fi

TAU=$(echo "scale=6;($SQRT_TAU)^2" | bc)
ONE_OVER_SQRT_TAU=$(echo "scale=6;1./$SQRT_TAU" | bc)




# Make a list of physics model parameters and set output names for the T2W rootfile and the postfit rootfile

if [ "$do_recoFit" = true ] ; then
    PHYSICSMODELPARAMETERS=$( python $CMSSW_BASE/src/scripts/getPhysicsModelParameterStrings.py $NBINS $NCATS )
    DCroot="${DC%.txt}_recoMuFit.root"
    DCrootpostfit="${DC%.txt}_recoMuFit_postfit.root"
else
    PHYSICSMODELPARAMETERS=$( python $CMSSW_BASE/src/scripts/getPhysicsModelParameterStrings.py $NBINS )
    DCroot="${DC%.txt}_genMuFit.root"
    DCrootpostfit="${DC%.txt}_genMuFit_postfit.root"
fi


########################################
# Run
########################################

if [ "$do_pT" = true ] ; then    
    cd pt_moriond17
else
    cd nJets_moriond17
fi


if [ "$do_T2WS" = true ] ; then

    chapter "Running text2workspace"
    echo "Datacard: $DC"
    echo "Output:   $DCroot"
    echo


    if [ "$do_pT" = true ] ; then    

        if [ "$do_recoFit" = false ] ; then
            # pT genFit
            text2workspace.py \
                $DC \
                -o $DCroot \
                -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
                --PO verbose \
                --PO 'map=.*/InsideAcceptance_genPt_0p0to15p0:r0[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_15p0to30p0:r1[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_30p0to45p0:r2[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_45p0to85p0:r3[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_85p0to125p0:r4[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_125p0to200p0:r5[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genPt_200p0to10000p0:r6[1,0,2]' \

        elif [ "$do_recoFit" = true ] ; then
            # pT recoFit
            text2workspace.py \
                $DC \
                -o $DCroot \
                -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
                --PO verbose \
                --PO 'map=.*SigmaMpTTag_0_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat2[1,0,2]'
        fi


    else

        if [ "$do_recoFit" = false ] ; then
            # nJets genFit
            text2workspace.py \
                $DC \
                -o $DCroot \
                -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
                --PO verbose \
                --PO 'map=.*/InsideAcceptance_genNjets2p5_m0p5to0p5:r0[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genNjets2p5_0p5to1p5:r1[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genNjets2p5_1p5to2p5:r2[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genNjets2p5_2p5to3p5:r3[1,0,2]' \
                --PO 'map=.*/InsideAcceptance_genNjets2p5_3p5to100p0:r4[1,0,2]'


        elif [ "$do_recoFit" = true ] ; then
            # nJets genFit
            text2workspace.py \
                $DC \
                -o $DCroot \
                -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
                --PO verbose \
                --PO 'map=.*SigmaMpTTag_0_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_0_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat0[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_1_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat1[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat2[1,0,2]' \
                --PO 'map=.*SigmaMpTTag_2_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat2[1,0,2]'

        fi
    fi
fi



if [ "$do_bestfit" = true ] ; then

    chapter "Running the best fit (asimov)"
    echo "Datacard: $DCroot"
    echo "Output:   $DCrootpostfit"
    echo

    combine \
        $DCroot \
        -M MultiDimFit \
        -t \
        -1 \
        --saveWorkspace \
        -n asimov_fit \
        --setPhysicsModelParameters $PHYSICSMODELPARAMETERS,OneOverSqrtTau=$ONE_OVER_SQRT_TAU \
        -m 125 \
        --minimizerStrategy 2
    mv higgsCombineasimov_fit.MultiDimFit.mH125.root $DCrootpostfit

fi


if [ "$do_generateToys" = true ] ; then

    chapter "Generating toys"
    echo "Asimov fit: $DCrootpostfit"
    echo

    allBkgVariables=$(python ../getCommaSeparatedListOfVariablesToFreeze.py $DCrootpostfit)

    combine \
        $DCrootpostfit \
        -M GenerateOnly \
        --setPhysicsModelParameters $PHYSICSMODELPARAMETERS,OneOverSqrtTau=$ONE_OVER_SQRT_TAU \
        -n asimov_toy \
        -m 125  \
        --saveToys \
        -t $NTOYS \
        --toysFrequentist --bypassFrequentistFit \
        --freezeNuisances $allBkgVariables \
        --snapshotName MultiDimFit

fi


if [ "$do_covMat" = true ] ; then

    allBkgVariables=$(python ../getCommaSeparatedListOfVariablesToFreeze.py $DCrootpostfit)

    combine \
        $DCrootpostfit \
        -M MultiDimFit \
        -n covmat_asimov_toy \
        -m 125 \
        -t $NTOYS \
        --toysFrequentist --bypassFrequentistFit \
        --algo none \
        --snapshotName MultiDimFit \
        --setPhysicsModelParameters $PHYSICSMODELPARAMETERS,OneOverSqrtTau=$ONE_OVER_SQRT_TAU \
        --freezeNuisances $allBkgVariables \
        -v 2
fi



cd -


