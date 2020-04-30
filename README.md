# ENA-metadata
R code to analyze the metadata downloaded from the ENA advanced search. Preparation of manuscript about data mining and status of ENA metadata.

### Related documents and files
* [google doc](https://docs.google.com/document/d/17AgTy6IJeaA7_hHMbWa2i5jO3LQPrHuYP1UDupiRFic/edit?usp=sharing) for manuscript preparation
* [owncloud](https://zmtcloud.zmt-bremen.de/index.php/s/Kwu0ty1bHIoJamK) with input data and associated references (password: metadataENA)

## Search queries
### amplicon data sets on ENA
On 05.03.2020, the [ENA advanced search](https://www.ebi.ac.uk/ena/browser/advanced-search) was used to retrieve the metadata (all available fields) for **raw reads** with the following search query:
```
tax_tree(410657) AND instrument_platform="ILLUMINA" AND library_layout="PAIRED" AND library_selection="PCR" AND library_strategy="AMPLICON" AND library_source="METAGENOMIC"
```

This search was updated on 28.04.2020. As the planned short commentary about the status of metadata on ENA is independent of the data mining studies, which were the inspiration for this paper, the results of the updated search will be used.

Search strategy:
```
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=instrument_platform%3D%22ILLUMINA%22%20AND%20library_layout%3D%22PAIRED%22%20AND%20library_selection%3D%22PCR%22%20AND%20library_strategy%3D%22AMPLICON%22%20AND%20library_source%3D%22METAGENOMIC%22&fields=accession%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbroker_name%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccompleteness_score%2Ccontamination_score%2Ccountry%2Ccram_index_aspera%2Ccram_index_ftp%2Ccram_index_galaxy%2Ccultivar%2Cculture_collection%2Cdepth%2Cdescription%2Cdev_stage%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_package%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_title%2Cexperimental_factor%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cgermline%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_layout%2Clibrary_name%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocation%2Clon%2Cmating_type%2Cnominal_length%2Cnominal_sdev%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Crun_accession%2Crun_alias%2Csalinity%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_sex%2Ctarget_gene%2Ctax_id%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Cvariety&limit=1000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"
```

Later, these results will be filtered by taxid, either to calculate separate statistics or to focus only on ecological metagenomes, which will be the main target for the analysis.


### WGS data sets on ENA
