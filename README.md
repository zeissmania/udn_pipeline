### Introduction

This package is designed for automating the UDN analysis process.

The scripts are written in python3



#### dependency

No package dependency.

The only dependency is the annotSV (if you need to run the annotatioin using this pipeline)

you could either specify the local annotSV path by

```
--annotsv  /path/to/annotSV
```

or run this task on ACCRE, the annotSV docker file is already on ACCRE, you don't need to specify the path.



#### 

### Files in need

1. VCF file, 

   the filename should be like   

   ```
   *UDN123456[_-][PMSDF]*vcf.gz
   ```

   for example, `919160-UDN217312-D.slm.cnv.vcf.gz`

2. phenotype file 

    filename  should ends with   "_terms.txt"  

   one phenotype per line



### Usage

1. run the whole pipeline

```
python3  run_udn.py   <project_name>  [other options]
```

2. only  get the AMELIE information

```
python amelie_api.py <project_name>  <file_gene_list>  <file_pure_hpo_id>
```



### File output

for example, the project_name is "case_test"

#### amele_api.py

1. `case_test.amelie.lite.txt`

   the AMELIE score for each gene

2. `case_test.amelie.parsed.pdict`  and `case_test.amelie.pdict` 

   the json information got from AMELIE API

3. `case_test.amelie_final_result.full.txt`

   the  detailed information for each gene, ( paper and related information)

#### run_udn.py

1. `case_test_terms_pure_hpo.txt`   and  `case_test_terms_hpo.txt`

   convert phenotype terms to HPO ID

2. `case_test_P_annotated.txt`   

   the original annotation file got from annotSV

   `_P_`  could be also  M/F/D/S  ( s for sibling,  M for mother, F/D for father)

   

3. `case_test_P_annotated.filtered.txt`

   basic filtering for the annotated result

   a) FILTER == PASS

   b)  only for proband annotation result:   gnomad MAF < 0.01,  annotAV ranking >=2

   

   `_P_`  could be also  M/F/D/S  ( s for sibling,  M for mother, F/D for father)

4. `case_test_P.extracted_fields.txt`

   extracted useful fields from `case_test_P_annotated.filtered.txt`

   `_P_`  could be also  M/F/D/S  ( s for sibling,  M for mother, F/D for father)

5. `case_test_D.filtered.genelist`

   the gene list after filtering, used for AMELIE 

6. `project_name_final_result.sorted.tsv`

   the main result, include  all the information for the final report

   need to select manually based on the `case_test.amelie_final_result.full.txt`  and  OMIM 

   genes sorted by AMELIE score

   

