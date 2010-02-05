k=0
folder="CosmicsBeamCommissioning09-PromptReco-v1RECO"
run=`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/StreamExpressCRAFT09-RpcCalHLT-v1ALCARECO/ | grep 1 -c`
rm list.php
wget http://higgs.skku.ac.kr/test/list.php

for j in `seq 1 ${run}`
do

   runnumber=$((`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/ | sed -n -e $((${j}))p | grep 1`))
   
   is=not 
   web=`grep -c 1 list.php`
   for l in `seq 1 $((${web}+1))`
   do
       web1=$((`cat list.php | sed -n -e $((${l}))p | grep 1`))
       if [ $runnumber = $web1 ]; then
          is=is
       fi
   done
   if [ $is = "not" ]; then

      gif=`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/ | grep -c gif`
      if [ "1" = "${gif}" ]; then
         k=$((${k}+1))
         Ngif=$((`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/$folder/${runnumber}/gif/ | grep -c gif`))
         echo 'Run: '${runnumber}', : '${Ngif}' : '$k
         zero=$((0))
	 if [ $Ngif -gt $zero ]; then
              echo "copying run ${runnumber} from castor, submitting to server..."
             ./getrun.sh ${runnumber}
         fi

                   
       #  if[ "0" == "${Ngif}" ]; then

#            mkdir run${runnumber}

#            for i in `seq 1 ${Ngif}`
#            do
#               rfcp /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/StreamExpressCRAFT09-RpcCalHLT-v1ALCARECO/${runnumber}/gif/`nsls /castor/cern.ch/cms/store/caf/user/ccmuon/RPC/GlobalRuns/StreamExpressCRAFT09-RpcCalHLT-v1ALCARECO/${runnumber}/gif/ | sed -n -e $((${i}))p`  ./run${runnumber}/
#            done
    #     fi
      fi
   fi
done
