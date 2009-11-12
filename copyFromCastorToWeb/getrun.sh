
   runnumber=${1}
#   folder="CosmicsBeamCommissioning09-PromptReco-v1RECO"
folder="BeamHaloBeamCommissioning09-PromptReco-v1RECO"
#   gif=`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/StreamExpressCRAFT09-RpcCalHLT-v1ALCARECO/${runnumber}/ | grep -c gif`
   gif=`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/ | grep -c gif`

   if [ "1" == "${gif}" ]; then
      Ngif=`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/gif/ | grep -c gif`
      echo 'Run: '${runnumber}', : '${Ngif}

    #  if[ "0" == "${Ngif}" ]; then

         mkdir run${runnumber}
	./mkdirAtServer.py ${runnumber}

         for i in `seq 1 ${Ngif}`
         do
           rfcp /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/gif/`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/gif/ | sed -n -e $((${i}))p`  ./run${runnumber}/
           ./submitToWeb.py $runnumber `pwd`/run${runnumber}/`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/gif/ | sed -n -e $((${i}))p`
         done
    #  fi
   fi

