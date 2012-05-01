
# start by getting the descendants of the Vira ancestral node
grep -w "10239" data/nodes.dmp  | awk '{print $1}' | uniq > tmp.virus.txt

SIZE_OLD=`wc -l tmp.virus.txt | awk '{print $1}'`
SIZE_NEW=0

echo $SIZE_OLD
ITER=0
while [ $SIZE_OLD != $SIZE_NEW ]
do
	ITER=$(expr $ITER + 1)
	DATE=`date`
	echo $DATE interation [ $ITER ]
	SIZE_OLD=`wc -l tmp.virus.txt | awk '{print $1}'`
	echo      -> old size [ $SIZE_OLD ]
	grep -wf tmp.virus.txt data/nodes.dmp  | awk '{print $1}'  >> tmp.virus.txt
	uniq tmp.virus.txt > tmp.virus2.txt
	mv tmp.virus2.txt tmp.virus.txt
	SIZE_NEW=`wc -l tmp.virus.txt | awk '{print $1}'`
	echo      -> new size [ $SIZE_NEW ]
done

