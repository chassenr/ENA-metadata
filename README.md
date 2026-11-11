# ENA-metadata

Assessment of metadata quality for metabarcoding sequencing data sets from environmental samples.

### Search strategy

**Curl request:**

```
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(12908)%20AND%20(instrument_platform%3D%22illumina%22%20AND%20library_layout%3D%22paired%22%20AND%20library_strategy%3D%22AMPLICON%22%20AND%20library_selection%3D%22PCR%22%20AND%20(%20library_source%3D%22genomic%22%20OR%20library_source%3D%22metagenomic%22%20))&fields=run_accession%2Cexperiment_title%2Ctax_id%2Cinstrument_platform%2Cage%2Caligned%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbam_aspera%2Cbam_bytes%2Cbam_ftp%2Cbam_galaxy%2Cbam_md5%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbisulfite_protocol%2Cbroad_scale_environmental_context%2Cbroker_name%2Ccage_protocol%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Cchip_ab_provider%2Cchip_protocol%2Cchip_target%2Ccollected_by%2Ccollection_date%2Ccollection_date_end%2Ccollection_date_start%2Ccompleteness_score%2Ccontamination_score%2Ccontrol_experiment%2Ccountry%2Ccultivar%2Cculture_collection%2Cdatahub%2Cdepth%2Cdescription%2Cdev_stage%2Cdisease%2Cdnase_protocol%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_medium%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_target%2Cexperimental_factor%2Cexperimental_protocol%2Cextraction_protocol%2Cfaang_library_selection%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfile_location%2Cfirst_created%2Cfirst_public%2Cgermline%2Chi_c_protocol%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_scientific_name%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_construction_protocol%2Clibrary_gen_protocol%2Clibrary_layout%2Clibrary_max_fragment_size%2Clibrary_min_fragment_size%2Clibrary_name%2Clibrary_pcr_isolation_protocol%2Clibrary_prep_date%2Clibrary_prep_date_format%2Clibrary_prep_latitude%2Clibrary_prep_location%2Clibrary_prep_longitude%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocal_environmental_context%2Clocation%2Clocation_end%2Clocation_start%2Clon%2Cmarine_region%2Cmating_type%2Cncbi_reporting_standard%2Cnominal_length%2Cnominal_sdev%2Cpcr_isolation_protocol%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Cread_strand%2Crestriction_enzyme%2Crestriction_enzyme_target_sequence%2Crestriction_site%2Crna_integrity_num%2Crna_prep_3_protocol%2Crna_prep_5_protocol%2Crna_purity_230_ratio%2Crna_purity_280_ratio%2Crt_prep_protocol%2Crun_alias%2Crun_date%2Csalinity%2Csample_accession%2Csample_alias%2Csample_capture_status%2Csample_collection%2Csample_description%2Csample_material%2Csample_prep_interval%2Csample_prep_interval_units%2Csample_storage%2Csample_storage_processing%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_project%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_date%2Csequencing_date_format%2Csequencing_location%2Csequencing_longitude%2Csequencing_method%2Csequencing_primer_catalog%2Csequencing_primer_lot%2Csequencing_primer_provider%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstatus%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmission_tool%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_read_type%2Ctag%2Ctarget_gene%2Ctax_lineage%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Ctransposase_protocol%2Cvariety&limit=10000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" > amplicon_results_read_run_tsv.txt
```

**Search criteria:**

```
tax_tree(12908) AND (instrument_platform="illumina" AND library_layout="paired" AND library_strategy="AMPLICON" AND library_selection="PCR" AND ( library_source="genomic" OR library_source="metagenomic" ))
```

The search was conducted on 2025-09-02. The search query included data sets from all unclassified sequence entries, since amplicon data sets from environmental samples may have an NCBI taxon ID other than 'metagenomes' (and children). Futhermore, we included both METAGENOMIC and GENOMIC as library source to account for potentially wrongly submitted data. This was done to quantify not just incomplete but also inconsistent metadata. This search strategy was based on personal experience in manual data curation during data mining studies and conversations about data archiving practices within the scientific community.


### S2ORC download

* Check the date of the newest release here: `https://api.semanticscholar.org/datasets/v1/release/latest`
* Get a new API key (it expires after ca. 60 days)
* Use [download_s2orc.py](https://github.com/chassenr/ENA-metadata/blob/master/download_s2orc.py) to download all the files (replace API key in line 8). 
* The required dependencies should be included in the base environment of Anaconda, but please see [here](https://github.com/chassenr/ENA-metadata/blob/master/base_env.yaml).

```
RELEASE_DATE = `2025-08-19`  
API_BASE = `https://api.semanticscholar.org/datasets/v1/release/2025-08-19/dataset/s2orc`
```

### Download full PMC dump

* Create a directory for the PMC download: ```PMC_data```
* Enter this directory and download all Open Access PMC articles: 

```
lftp -c "open ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/; mirror --parallel=4 ."
lftp -c "open ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/xml/; mirror --parallel=4 ."
lftp -c "open ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/xml/; mirror --parallel=4 ."
```

:warning: NCBI changed access to PMC from ftp to AWS in 2026. The above documentation still includes the ftp download instructions since the PMC corpus was accessed on 2025-09-02. To repeat the analysis, see [here] for updated instructions.


### File structure

The main directory of the analysis on the IOW servers is: ```/bio/Analysis_data/lost_found/tmp_hassenru/ENA_data_mining```. This directory contains several subdirectories:

```
ENA_search_tsv   # TSV output of the ENA advanced search
ENA_search_xml   # Study, sample and experiment XML of all entries in the search output
ENA_xref_out     # XREF output listing publications linked to studies in search output
ENVO_files       # Files (and parsed information) from the Environmental Ontology
PMC_data         # PMC data download
R_analysis       # Intermediate files produced during the analysis in R
s2orc_data       # S2ORC data download
Summaries        # Summary tables used for figure generation
```


### Filtering initial metadata to remove non-target hits

For the assessment of metadata quality for environmental metabarcoding sequencing data on ENA, human material should be excluded as well as samples from clinical settings. This filtering step relies only on the NCBI taxon ID as a mandatory field using a controlled vocabulary. The filter was implemented using this interactive [R script](https://github.com/chassenr/ENA-metadata/blob/master/status_ena_illumina_amplicon.R). The R environment used for the analysis is defined [here](https://github.com/chassenr/ENA-metadata/blob/master/r-4.0.0_env.yaml)


### File download and processing

The download of the sample XML was unsuccessful from the ENA advanced search results website. Downloading individual xml independently:

```
cd ENA_search_xml
mkdir sample_xml
cut -f84 ../ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt | sed '1d' | sort -u | parallel -j32 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O sample_xml/{}.xml'
# repeat for those where download failed (404 error)
cd sample_xml
find . -empty > ../empty_sample.txt
cat ../empty_sample.txt | xargs -n1 basename | sed 's/\.xml//' | parallel -j10 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O {}.xml'
find . -empty > ../empty_check.txt
diff ../empty_sample.txt ../empty_check.txt
mv ../empty_check.txt ../empty_sample.txt
# instead of using find, ls -l + awk to identify empty files is so much faster
# there are 206 samples for which information could now be retrieved...
```
I am not checking biosamples for the missing information considering that only a minor fraction of samples are affected.

Sample xml files were converted to tab-separated text for further processing using [this script](https://github.com/chassenr/ENA-metadata/blob/master/convert_xml_sample_attributes_to_csv.py).

```
cd ENA_search_xml/sample_xml
mkdir ../sample_csv
conda activate base
# replace /path/to/repo with file location of this repository
ls -1 | parallel -j96 'python /path/to/repo/convert_xml_sample_attributes_to_csv.py {} ../sample_csv'
```

Parsed csv files were combined into 1 table by bind_cols in [R](https://github.com/chassenr/ENA-metadata/blob/master/status_ena_illumina_amplicon.R) and column headers were manually screened based on their name and occasional cross-checks of their values whether they may contain information relevant for specific metadata.


Retrieving study-level information from ENA:
```
cd ENA_search_xml
mkdir study_xml
cut -f108 ../ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt | sed '1d' | sort -u | parallel -j32 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O study_xml/{}.xml'
# repeat for those where download failed
cd study_xml
find . -empty > ../empty_study.txt
cat ../empty_study.txt | xargs -n1 basename | sed 's/\.xml//' | parallel -j32 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O {}.xml'
find . -empty > ../empty_check.txt
diff ../empty_study.txt ../empty_check.txt
mv ../empty_check.txt ../empty_study.txt
# repeat as needed
```
There are some (66) studies whose accession numbers do not exist on ENA.


Retrieve publication via XREF if available this way:
```
# xref on study level
cd ENA_xref_out
mkdir xref_primary
cut -f108 ../ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt | sed '1d' | sort -u | parallel -j10 'wget https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession={} -O xref_primary/{}.txt > xref_primary/{}.log 2>&1'
grep "HTTP request sent" xref_primary/*.log | grep -c -v "200" # all good
cat xref_primary/*.txt > xref_results_primary_study.txt
cut -f1 xref_results_primary_study.txt | sort | uniq -c
grep -E "^EuropePMC|^PubMed" xref_results_primary_study.txt > xref_results_primary_study_pmc.txt
# retrieve additional info from secondary study accession number
mkdir xref_secondary
cut -f95 ../ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt | sed '1d' | sort -u | parallel -j1 'wget https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession={} -O xref_secondary/{}.txt > xref_secondary/{}.log 2>&1'
grep "HTTP request sent" xref_secondary/*.log | grep -c -v "200" # all good
cat xref_secondary/*.txt > xref_results_secondary_study.txt
cut -f1 xref_results_secondary_study.txt | sort | uniq -c
grep -E "^EuropePMC|^PubMed" xref_results_secondary_study.txt > xref_results_secondary_study_pmc.txt
# check redundancy between primary and secondary accession XREF hits
cut -f6,7 xref_results_primary_study_pmc.txt | tr '\t' '\n' | sort -u | sed '/^$/d' | grep -c -w -v -F -f - xref_results_secondary_study_pmc.txt
# there are no hits based on the secondary accession that are not also part of the primary accession output
```

Since XREF may not find all data set citations, search PMC and S2ORC for regular expressions corresponding to primary and secondary study accessions.
```
# pmc
cd R_analysis
grep -HEorw "PRJ[A-Z]{2}[0-9]{4,6}" ../PMC_data | awk -F/ '{print $NF}' > pmc_study_pid_numbers.txt
grep -HEorw "[A-Z]{1}RP[0-9]{6,8}" ../PMC_data | awk -F/ '{print $NF}' > pmc_study_sid_numbers.txt
cat pmc_study_* | cut -d':' -f2 | sort -u > pmc_matches_study_accession.txt
# s2orc
# run on phy-3 since jq was available system-wide (replace /path/to/repo with file location of this repository)
ls -1 ../s2orc_data/ | parallel -j96 '/path/to/repo/s2orc_bash_script_pid.sh ../s2orc_data/{} s2orc_pid_out/{}_pid.txt'
ls -1 ../s2orc_data/ | parallel -j96 '/path/to/repo/s2orc_bash_script_sid.sh ../s2orc_data/{} s2orc_sid_out/{}_sid.txt'
cat s2orc_*/*.txt | cut -d' ' -f2 | sort -u > s2orc_matches_study_accession.txt
```


There may be further information about the amplified region in the experiment XML in addition to the library construction protocol (which is provided in the tsv output) under DESIGN_DESCRIPTION and EXPERIMENT_ATTRIBUTES.

```
cd ENA_search_xml
mkdir experiment_xml
cut -f37 ../ENA_search_tsv/amplicon_results_read_run_tsv_non_human.txt | sed '1d' | sort -u | parallel -j20 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O experiment_xml/{}.xml'
# repeat for those where download failed (404 error)
cd experiment_xml
ls -l > ../list_exp.txt
sed -E 's/ +/\t/g' ../list_exp.txt | awk '$5 == 0' | cut -f9 | sed 's/\.xml//' | parallel -j20 'wget https://www.ebi.ac.uk/ena/browser/api/xml/{} -O {}.xml'
```

Experiment xml files were converted to tab-separated text for further processing using [this script](https://github.com/chassenr/ENA-metadata/blob/master/convert_xml_experiment_attributes_to_csv.py).

```
cd ENA_search_xml/experiment_xml
mkdir ../experiment_csv
conda activate base
# replace /path/to/repo with file location of this repository
ls -1 | parallel -j48 'python /path/to/repo/convert_xml_experiment_attributes_to_csv.py {} ../experiment_csv'
```


### Metadata quality assessment

For each of the selected metadata fields, the code for their assessment is provided in this [Rscript](https://github.com/chassenr/ENA-metadata/blob/master/status_ena_illumina_amplicon.R), which also includes the chi-square statistics.


### Figure generation

Data visualization is described in [figure_generation.R](https://github.com/chassenr/ENA-metadata/blob/master/figure_generation.R).