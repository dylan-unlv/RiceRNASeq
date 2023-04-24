import json
import os
import pandas as pd

colnames = ['Go_Term','label','module','FDR','p-value','fold-enrichment','plus-minus']
agg_data = pd.DataFrame(columns = colnames)


for filename in os.listdir('data/panther_out'):
    jdata = json.load(open('data/panther_out/'+filename,'r'))
    for result in jdata['results']['result']:
        if (result['fdr'] < 0.05) and ('id' in result['term']):
            row = pd.DataFrame([[result['term']['id'], result['term']['label'], filename.split('.')[0],result['fdr'], result['pValue'],result['fold_enrichment'],result['plus_minus']]], columns=colnames)
            agg_data = agg_data.append(row)

agg_data.sort_values(by = ['FDR'], axis = 0)
agg_data.to_csv('data/aggregated.tsv',sep='\t',index=False)

