# CURL request:
# curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(410657)%20AND%20library_selection%3D%22PCR%22%20AND%20library_strategy%3D%22AMPLICON%22%20AND%20library_layout%3D%22PAIRED%22%20AND%20instrument_platform%3D%22ILLUMINA%22%20AND%20library_source%3D%22METAGENOMIC%22&fields=accession%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbroker_name%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccompleteness_score%2Ccontamination_score%2Ccountry%2Ccram_index_aspera%2Ccram_index_ftp%2Ccram_index_galaxy%2Ccultivar%2Cculture_collection%2Cdepth%2Cdescription%2Cdev_stage%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_package%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_title%2Cexperimental_factor%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cgermline%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_layout%2Clibrary_name%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocation%2Clon%2Cmating_type%2Cnominal_length%2Cnominal_sdev%2Cparent_study%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Crun_accession%2Crun_alias%2Csalinity%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_sex%2Ctarget_gene%2Ctax_id%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Cvariety&limit=500000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"

# sample xml can already be downloaded from search website

# get xml for experiment and run
# work on olorin and 'borrow' Mara's taxdb env for kraken2 (which contains parallel)
conda activate  /storage/hdd2/mmaeke/Github/taxdb-integration/.snakemake/conda/bcda7581
cut -f1,33,77 results_read_run_tsv.txt > results_read_run_tsv_sample_exp_run.txt
sed '1d' results_read_run_tsv_sample_exp_run.txt | cut -f2 | sed 's/^/https\:\/\/www\.ebi\.ac\.uk\/ena\/browser\/api\/xml\//' | parallel -j200 'wget -qO- {}' | sed '/<\/EXPERIMENT_SET>/,/<EXPERIMENT_SET>/d' >> experiment.xml
# some entries are suppressed (no idea why)
# append </EXPERIMENT_SET>
cp experiment.xml experiment_ori.xml
echo "</EXPERIMENT_SET>" >> experiment.xml

sed '1d' results_read_run_tsv_sample_exp_run.txt | cut -f3 | sed 's/^/https\:\/\/www\.ebi\.ac\.uk\/ena\/browser\/api\/xml\//' | parallel -j200 'wget -qO- {}' | sed '/<\/RUN_SET>/,/<RUN_SET>/d' >> run.xml
# same entries are suppressed (no idea why)
cp run.xml run_ori.xml
echo "</RUN_SET>" >> run.xml

# filter ENA runs based on correct ontology terms for biome, material, feature
# check for biome ID
cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > biome_id.txt
# check for biome name
cut -f1 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > biome_name.txt
# check for biome name (relaxed)
cut -f4 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > biome_name2.txt

# check for material ID
cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > material_id.txt
# check for material name
cut -f2 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > material_name.txt
# check for material name (relaxed)
cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > material_name2.txt

# check for feature ID
cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f3 envo_formatted.txt) > feature_id.txt
# check for feature name
cut -f3 ena_out_envo.txt | grep -n -a -f <(cut -f5 envo_formatted.txt) > feature_name.txt
# check for feature name (relaxed)
cut -f4 ena_out_envo.txt | grep -n -a -f <(cut -f6 envo_formatted.txt) > feature_name2.txt


# Retrieve Pubmed ID via XREF
# run on server:
mkdir ENA_xref
cd ENA_xref

# check with primary study accession
mkdir xref_primary_study_accnos
cd xref_primary_study_accnos
cut -f1 ../../ena_subset_sulfurimonas.txt | sed '1d' | while read line
do
  curl -s -o ${line}"_xref.txt" "https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession="${line}"&source=EuropePMC"
done

# filter studies with no hits (file size is 17)
ls -l *xref.txt | sed 's/ \+/\t/g' | cut -f5 | sort | uniq -c
find . -maxdepth 1 -size -20c | grep "xref" > ../no_hits_primary.txt
# 5065 no hits
find . -maxdepth 1 -size -20c | grep "xref" | xargs rm
cd ..

# get secondary accnos for these studies
# curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=study&fields=breed%2Cbroker_name%2Ccenter_name%2Ccultivar%2Cfirst_public%2Cgeo_accession%2Cisolate%2Ckeywords%2Clast_updated%2Cparent_study%2Cscientific_name%2Csecondary_study_accession%2Cstrain%2Cstudy_accession%2Cstudy_description%2Cstudy_name%2Cstudy_title%2Ctax_id&includeAccessionType=study&includeAccessions=PRJEB10576%2CPRJEB11384%2CPRJEB14127%2CPRJEB18774%2CPRJEB31776%2CPRJEB7448%2CPRJNA274359%2CPRJNA282077%2CPRJNA299110%2CPRJNA318932%2CPRJNA330786%2CPRJNA331054%2CPRJNA341261%2CPRJNA349764%2CPRJNA360358%2CPRJNA414441%2CPRJNA434752%2CPRJNA485064%2CPRJNA498402%2CPRJNA524261%2CPRJNA549457%2CPRJNA563517%2CPRJNA564579&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"

# check with secondary study accession
mkdir xref_secondary_study_accnos
cd xref_secondary_study_accnos
cut -f12 ../study_search_results_no_hits_primary.txt | sed '1d' | while read line
do
  curl -s -o ${line}"_xref.txt" "https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession="${line}"&source=EuropePMC"
done

# filter studies with no hits (file size is 17)
find . -maxdepth 1 -size -20c | grep "xref" > ../no_hits_secondary.txt
# 4729 no hits
find . -maxdepth 1 -size -20c | grep "xref" | xargs rm
cd ..

# get primary for no XREF hits
sed -e 's/\.\///' -e 's/_xref\.txt//' no_hits_secondary.txt | grep -F -f - study_search_results_no_hits_primary.txt | cut -f14 > no_xref_study_accnos.txt

