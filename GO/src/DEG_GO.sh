#!/bin/bash

#start with upregulated
DEFILE=$(ls data/DEG/upregulated/)  


for F in $DEFILE
do
	IDs=$(cut -f1 'data/DEG/upregulated/'$F | tail -n +2 | head )
	treatment=$(echo $F | cut -f2 -d'_' | cut -f1 -d'.')	
	#echo $treatment
	#echo $IDs
	curl -X POST "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="${IDs}"&organism=39947&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json" > data/DEG/panther_out/upregulated/${treatment}.BP.json
	curl -X POST "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="${IDs}"&organism=39947&annotDataSet=GO%3A0003674&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json" > data/DEG/panther_out/upregulated/${treatment}.up.MF.json
	curl -X POST "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="${IDs}"&organism=39947&annotDataSet=GO%3A0005575&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json" > data/DEG/panther_out/upregulated/${treatment}.up.CC.json

done
	
