# ENA-metadata
R code to analyze the metadata downloaded from the ENA advanced search. Preparation of manuscript about data mining and status of ENA metadata.

## Search queries

### amplicon data sets on ENA
On 13.12.2020, the [ENA advanced search](https://www.ebi.ac.uk/ena/browser/advanced-search) was used to retrieve the metadata (all available fields) for **raw reads** with the following search query:

```
tax_tree(410657) AND library_selection = "PCR" AND library_strategy = "AMPLICON" AND library_layout = "PAIRED" AND instrument_platform = "ILLUMINA" AND library_source = "METAGENOMIC"
```

Search strategy:

```
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(410657)%20AND%20library_selection%3D%22PCR%22%20AND%20library_strategy%3D%22AMPLICON%22%20AND%20library_layout%3D%22PAIRED%22%20AND%20instrument_platform%3D%22ILLUMINA%22%20AND%20library_source%3D%22METAGENOMIC%22&fields=accession%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbroker_name%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccompleteness_score%2Ccontamination_score%2Ccountry%2Ccram_index_aspera%2Ccram_index_ftp%2Ccram_index_galaxy%2Ccultivar%2Cculture_collection%2Cdepth%2Cdescription%2Cdev_stage%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_package%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_title%2Cexperimental_factor%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cgermline%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_layout%2Clibrary_name%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocation%2Clon%2Cmating_type%2Cnominal_length%2Cnominal_sdev%2Cparent_study%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Crun_accession%2Crun_alias%2Csalinity%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_sex%2Ctarget_gene%2Ctax_id%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Cvariety&limit=500000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"
```


### WGS data sets on ENA
On 08.01.2021, the [ENA advanced search](https://www.ebi.ac.uk/ena/browser/advanced-search) was used to retrieve the metadata (all available fields) for **raw reads** with the following search query:

```
tax_tree(410657) AND library_selection = "RANDOM" AND library_strategy = "WGS" AND library_layout = "PAIRED" AND instrument_platform = "ILLUMINA" AND library_source = "METAGENOMIC" AND first_created<=2020-12-31
```

Only studies until the end of 2020 were considered.

Search strategy:

```
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(410657)%20AND%20instrument_platform%3D%22ILLUMINA%22%20AND%20library_layout%3D%22PAIRED%22%20AND%20library_selection%3D%22RANDOM%22%20AND%20library_source%3D%22METAGENOMIC%22%20AND%20library_strategy%3D%22WGS%22%20AND%20first_created%3C%3D2020-12-31&fields=accession%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbroker_name%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccompleteness_score%2Ccontamination_score%2Ccountry%2Ccram_index_aspera%2Ccram_index_ftp%2Ccram_index_galaxy%2Ccultivar%2Cculture_collection%2Cdepth%2Cdescription%2Cdev_stage%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_package%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_title%2Cexperimental_factor%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cgermline%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_layout%2Clibrary_name%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocation%2Clon%2Cmating_type%2Cnominal_length%2Cnominal_sdev%2Cparent_study%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Crun_accession%2Crun_alias%2Csalinity%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_sex%2Ctarget_gene%2Ctax_id%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Cvariety&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"
```

For both searches, results were downloaded in tsv and xml format for further curation and plotting.


## Analysis update

The above search was repeated on 12.07.2023. In the meantime the number of runs had tripled and the ENA API had been updated. The evaluation of the status of sequence metadata on ENA was repeated to account for those changes.
The download of the sample XML was unsuccessful from the ENA advanced search results website. Downloading individual xml independently:

```
# on bio-48
cd /home/hassenru/ENA_data_mining
cut -f2 results_read_run_tsv.txt | sed '1d' | sort | uniq | parallel -j32 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O sample_xml/{}.xml'
# for those where download failed, get sample metadata from biosamples
cd sample_xml
ls -l | sed '1d' | sed -r 's/ +/\t/g' | awk '$5 == 0' | cut -f9 | sed 's/\.xml//' | parallel -j32 ' wget https://www.ebi.ac.uk/biosamples/samples/{}.xml -O ../biosample_xml/{}.xml'
# this did not yield any useful results
# experiment xml
cd ..
cut -f46 results_read_run_tsv.txt | sed '1d' | sort | uniq | parallel -j32 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O experiment_xml/{}.xml'
```

Weirdly, there are many runs whose sample accession number does not exist (neither on samples nor on biosamples). Checking the biosample DB only retrieved additional records for 3 sample accessions. I am not including this output for now as it will not impact the overall results. Likewise there are also runs whose experiment accession does not exist (anymore).


Retrieve publication via XREF if available this way:

```
# xref on study level
cd xref_out
cut -f109 ../results_read_run_tsv.txt | sed '1d' | sort | uniq | while read line
do
  wget https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession=${line} -O tmp.txt
  sed '1d' tmp.txt >> xref_results_primary_study.txt
  rm tmp.txt
done
cut -f144 ../results_read_run_tsv.txt | sed '1d' | sort | uniq | while read line
do
  wget https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession=${line} -O tmp.txt
  sed '1d' tmp.txt >> xref_results_secondary_study.txt
  rm tmp.txt
done
```

Process parsed sample xml to extract all possible attribute names.

```
# on bio-48
cd /home/hassenru/ENA_data_mining/xml_parsed_csv/sample_csv
ls -1 | parallel -j30 "head -1 {} | tr '\t' '\n' " > ../all_sample_attribute_names.txt
cd ..
sort all_sample_attribute_names.txt | uniq -c | sed -e 's/^ *//' -e 's/ /\t/' > all_sample_attribute_names_counts.txt
```



