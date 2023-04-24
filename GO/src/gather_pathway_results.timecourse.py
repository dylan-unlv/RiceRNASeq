import json
import os
import pandas as pd

colnames = ['Go_Term','label','annotation','module','FDR','p-value','fold-enrichment','plus-minus','num_genes']
agg_data = pd.DataFrame(columns = colnames)

count = 0
for filename in os.listdir('data/DEG/results/upregulated_website/'):
    print(filename)
    if '.json' in filename:
        jdata = json.load(open('data/DEG/results/upregulated_website/'+filename,'r'))
        for result in jdata['overrepresentation']['group']:
            #i=0
            #rdict = result['result'][0]
            #print([rdict['term']['id'], rdict['term']['label'], filename.split('.')[0].split('_')[-1], result['result'][0]['input_list'][0]['list_name'].split('.')[0].split('_')[-1].strip('hr'), rdict['input_list'][i]['fdr'],rdict['input_list'][i]['pValue'] , rdict['input_list'][i]['fold_enrichment'], rdict['input_list'][i]['plus_minus'], rdict['input_list'][i]['number_in_list'])
            for rdict in result['result']:
                if not isinstance(rdict, str):
                    for i in range(0, len(rdict['input_list'])):
                        if (rdict['input_list'][i]['fdr'] < 0.05) and ('id' in rdict['term']):
                            row = pd.DataFrame([[rdict['term']['id'], rdict['term']['label'], filename.split('.')[0].split('_')[-1], result['result'][0]['input_list'][i]['list_name'].split('.')[0].split('_')[-1], rdict['input_list'][i]['fdr'],rdict['input_list'][i]['pValue'] , rdict['input_list'][i]['fold_enrichment'], rdict['input_list'][i]['plus_minus'], rdict['input_list'][i]['number_in_list']]], columns=colnames)
                            agg_data = agg_data.append(row)
                else:
                    for i in range(0, len(result['result']['input_list'])):
                        #print([[result['result']['term']['id'], result['result']['term']['label'], filename.split('.')[0].split('_')[1],result['result']['input_list'][i]['list_name'].split('.')[0].split('_')[-1].strip('hr'),result['result']['input_list'][i]['fdr'], result['result']['input_list'][i]['pValue'],result['result']['input_list'][i]['fold_enrichment'],result['result']['input_list'][i]['plus_minus'],result['result']['input_list'][i]['number_in_list']]]) 
                        if (result['result']['input_list'][i]['fdr'] < 0.05) and ('id' in result['result']['term']):
                            row = pd.DataFrame([[result['result']['term']['id'], result['result']['term']['label'], filename.split('.')[0].split('_')[1],result['result']['input_list'][i]['list_name'].split('.')[0].split('_')[-1],result['result']['input_list'][i]['fdr'], result['result']['input_list'][i]['pValue'],result['result']['input_list'][i]['fold_enrichment'],result['result']['input_list'][i]['plus_minus'],result['result']['input_list'][i]['number_in_list']]], columns=colnames)
                            agg_data = agg_data.append(row)
                    
                    

agg_data.sort_values(by = ['FDR'], axis = 0)
agg_data.to_csv('data/timecourse.up.aggregated.slim.tsv',sep='\t',index=False)

colnames = ['Go_Term','label','annotation','module','FDR','p-value','fold-enrichment','plus-minus','num_genes']
agg_data = pd.DataFrame(columns = colnames)


for filename in os.listdir('data/DEG/results/downregulated_website/'):
    print(filename)
    if '.json' in filename:
        jdata = json.load(open('data/DEG/results/downregulated_website/'+filename,'r'))
        for result in jdata['overrepresentation']['group']:
            #i=0
            #rdict = result['result'][0]
            #print([rdict['term']['id'], rdict['term']['label'], filename.split('.')[0].split('_')[-1], result['result'][0]['input_list'][0]['list_name'].split('.')[0].split('_')[-1].strip('hr'), rdict['input_list'][i]['fdr'],rdict['input_list'][i]['pValue'] , rdict['input_list'][i]['fold_enrichment'], rdict['input_list'][i]['plus_minus'], rdict['input_list'][i]['number_in_list'])
            for rdict in result['result']:
                if not isinstance(rdict, str):
                    for i in range(0, len(rdict['input_list'])):
                        if (rdict['input_list'][i]['fdr'] < 0.05) and ('id' in rdict['term']):
                            row = pd.DataFrame([[rdict['term']['id'], rdict['term']['label'], filename.split('.')[0].split('_')[-1], result['result'][0]['input_list'][i]['list_name'].split('.')[0].split('_')[-1], rdict['input_list'][i]['fdr'],rdict['input_list'][i]['pValue'] , rdict['input_list'][i]['fold_enrichment'], rdict['input_list'][i]['plus_minus'], rdict['input_list'][i]['number_in_list']]], columns=colnames)
                            agg_data = agg_data.append(row)
                else:
                    for i in range(0, len(result['result']['input_list'])):
                        #print([[result['result']['term']['id'], result['result']['term']['label'], filename.split('.')[0].split('_')[1],result['result']['input_list'][i]['list_name'].split('.')[0].split('_')[-1].strip('hr'),result['result']['input_list'][i]['fdr'], result['result']['input_list'][i]['pValue'],result['result']['input_list'][i]['fold_enrichment'],result['result']['input_list'][i]['plus_minus'],result['result']['input_list'][i]['number_in_list']]])
                        
                        if (result['result']['input_list'][i]['fdr'] < 0.05) and ('id' in result['result']['term']):
                            row = pd.DataFrame([[result['result']['term']['id'], result['result']['term']['label'], filename.split('.')[0].split('_')[1],result['result']['input_list'][i]['list_name'].split('.')[0].split('_')[-1],result['result']['input_list'][i]['fdr'], result['result']['input_list'][i]['pValue'],result['result']['input_list'][i]['fold_enrichment'],result['result']['input_list'][i]['plus_minus'],result['result']['input_list'][i]['number_in_list']]], columns=colnames)
                            agg_data = agg_data.append(row)
agg_data.sort_values(by = ['FDR'], axis = 0)
agg_data.to_csv('data/timecourse.down.aggregated.slim.tsv',sep='\t',index=False)
