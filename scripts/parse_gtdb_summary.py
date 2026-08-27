#!/usr/bin/env python3
import sys, json
for line in sys.stdin:
    d=json.loads(line)
    acc=d['accession']
    biosample=d['assembly_info']['biosample']['accession']
    sra_ids=[]
    for sid in d['assembly_info']['biosample'].get('sample_ids',[]):
        if sid.get('db')=='SRA':
            sra_ids.append(sid['value'])
    strain=d['assembly_info']['biosample'].get('strain','')
    print(acc, biosample, ','.join(sra_ids), strain, sep='\t')
