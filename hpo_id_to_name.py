
hpo1 = eval(open('/Users/files/ref/HPO/hpo_name_to_id_alt.pdict').read())
hpo2 = eval(open('/Users/files/ref/HPO/hpo_name_to_id.pdict').read())
hpo = {}
for k, v in hpo1.items():
    for hpoid in v:
        try:
            hpo[hpoid].append(k)
        except:
            hpo[hpoid] = [k]


for k, hpoid in hpo2.items():
    try:
        hpo[hpoid].append(k)
    except:
        hpo[hpoid] = [k]
