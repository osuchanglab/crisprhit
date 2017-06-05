#!/bin/bash

bmfile='benchmark.crisprhit'

if [[ ! -z "$1" ]]; then
	bmfile="$1"
fi

if [[ ! -e "$bmfile" ]]; then
	echo "Please provide an existing file. ($bmfile) does not exist."
fi

pattern='_D'

D_hq=`grep $pattern $bmfile | grep -c -P 'hq\t'`
D_hqp=`grep $pattern $bmfile | grep -c -P '\thq_priming'`
D_p=`grep $pattern $bmfile | grep -c -P '\tpriming'`
D_s=`grep $pattern $bmfile | grep -c -P '\tstable'`
D_o=`grep $pattern $bmfile | grep -c -P '\tother'`

pattern='_P'

P_hq=`grep $pattern $bmfile | grep -c -P 'hq\t'`
P_hqp=`grep $pattern $bmfile | grep -c -P '\thq_priming'`
P_p=`grep $pattern $bmfile | grep -c -P '\tpriming'`
P_s=`grep $pattern $bmfile | grep -c -P '\tstable'`
P_o=`grep $pattern $bmfile | grep -c -P '\tother'`

pattern='_S'

S_hq=`grep $pattern $bmfile | grep -c -P 'hq\t'`
S_hqp=`grep $pattern $bmfile | grep -c -P '\thq_priming'`
S_p=`grep $pattern $bmfile | grep -c -P '\tpriming'`
S_s=`grep $pattern $bmfile | grep -c -P '\tstable'`
S_o=`grep $pattern $bmfile | grep -c -P '\tother'`

echo -e "interaction\thq\thq_priming\tpriming\tstable\tother"
echo -e "DI\t$D_hq\t$D_hqp\t$D_p\t$D_s\t$D_o"
echo -e "Priming\t$P_hq\t$P_hqp\t$P_p\t$P_s\t$P_o"
echo -e "Stable\t$S_hq\t$S_hqp\t$S_p\t$S_s\t$S_o"
