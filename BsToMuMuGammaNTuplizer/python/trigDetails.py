import FWCore.ParameterSet.Config as cms

triggersOfInterest=cms.vstring( 
                                "HLT_DoubleMu4_3_Bs_v14",
                                "HLT_DoubleMu4_Jpsi_Displaced_v7",
                                "HLT_DoubleMu4_JpsiTrk_Displaced_v15",
                              )
bParkTriglist=cms.vstring( 
    "HLT_DoubleMu4_3_Bs",
    "HLT_DoubleMu4_Jpsi_Displaced",
    "HLT_DoubleMu4_JpsiTrk_Displaced_v15",
    'HLT_Mu10p5_IP3p5',
    'HLT_Mu12_IP6',
    'HLT_Mu7_IP4',
    'HLT_Mu8_IP3',
    'HLT_Mu8_IP5',
    'HLT_Mu8p5_IP3p5',
    'HLT_Mu9_IP5',
    'HLT_Mu9_IP6',
    )

bParkTriglist_full=cms.vstring( 
    'HLT_Mu12_IP6_ToCSCS_v1',
	'HLT_Mu12_IP6_part0_v2',
	'HLT_Mu12_IP6_part1_v2',
	'HLT_Mu12_IP6_part2_v2',
	'HLT_Mu12_IP6_part3_v2',
	'HLT_Mu12_IP6_part4_v2',
	'HLT_Mu9_IP5_ToCSCS_v1',
	'HLT_Mu9_IP5_part0_v2',
	'HLT_Mu9_IP5_part1_v2',
	'HLT_Mu9_IP5_part2_v2',
	'HLT_Mu9_IP5_part3_v2',
	'HLT_Mu9_IP5_part4_v2',
	'HLT_Mu7_IP4_ToCSCS_v1',
	'HLT_Mu7_IP4_part0_v2',
	'HLT_Mu7_IP4_part1_v2',
	'HLT_Mu7_IP4_part2_v2',
	'HLT_Mu7_IP4_part3_v2',
	'HLT_Mu7_IP4_part4_v2',
	'HLT_Mu9_IP4_ToCSCS_v1',
	'HLT_Mu9_IP4_part0_v2',
	'HLT_Mu9_IP4_part1_v2',
	'HLT_Mu9_IP4_part2_v2',
	'HLT_Mu9_IP4_part3_v2',
	'HLT_Mu9_IP4_part4_v2',
	'HLT_Mu8_IP5_ToCSCS_v1',
	'HLT_Mu8_IP5_part0_v2',
	'HLT_Mu8_IP5_part1_v2',
	'HLT_Mu8_IP5_part2_v2',
	'HLT_Mu8_IP5_part3_v2',
	'HLT_Mu8_IP5_part4_v2',
	'HLT_Mu8_IP6_ToCSCS_v1',
	'HLT_Mu8_IP6_part0_v2',
	'HLT_Mu8_IP6_part1_v2',
	'HLT_Mu8_IP6_part2_v2',
	'HLT_Mu8_IP6_part3_v2',
	'HLT_Mu8_IP6_part4_v2',
	'HLT_Mu9_IP6_ToCSCS_v1',
	'HLT_Mu9_IP6_part0_v3',
	'HLT_Mu9_IP6_part1_v3',
	'HLT_Mu9_IP6_part2_v3',
	'HLT_Mu9_IP6_part3_v3',
	'HLT_Mu9_IP6_part4_v3',
	'HLT_Mu8_IP3_ToCSCS_v1',
	'HLT_Mu8_IP3_part0_v3',
	'HLT_Mu8_IP3_part1_v3',
	'HLT_Mu8_IP3_part2_v3',
	'HLT_Mu8_IP3_part3_v3',
	'HLT_Mu8_IP3_part4_v3',
	'HLT_Mu9_IP0_part0_v2',
	'HLT_Mu9_IP3_part0_v2'
)

bParkTrigFilterModules=cms.vstring(
    'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP0Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP3Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q',
	'hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q'
   )