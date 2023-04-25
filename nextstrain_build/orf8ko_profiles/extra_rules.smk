ruleorder: filter_replace > filter

rule download_metadata_s3:
    message: "Downloading metadata from {params.address} -> {output.metadata}"
    output:
        metadata = "data/downloaded_meta.tsv.gz"
    params:
        address = "s3://nextstrain-ncov-private/metadata.tsv.gz"
    shell:
        """
        aws s3 cp {params.address} {output.metadata}
        """

rule add_metadata:
    message:
        """
        Deletion to metadata
        """
    input:
        metadata = rules.download_metadata_s3.output,
        #deletions = "data/washington_deletions_dedup.tsv"
        deletions = "data/washington_deletions.tsv"
    params:
        key = "strain",
        add = "Ns,gap,protein_length,orf8ko"
    output:
        combined = "data/gisaid_deletions.tsv.gz"
    shell:
        """
        gunzip -d -c {input.metadata} | tsv-join \
        --filter-file {input.deletions} \
         --key-fields {params.key} \
         --append-fields {params.add} \
         -H \
         --write-all ''\
         | gzip > {output.combined}
        """

rule include:
    message:
        """
        Adding all WA deletions to include
        """
    input:
        metadata = "results/{build_name}/metadata_with_index.tsv",
    output:
        include = "results/{build_name}/include.txt"
    shell:
        """
        tsv-filter -H --str-eq orf8ko:Yes {input.metadata} | cut -f 1 | sed -n '1!p' > {output.include}
        """

rule filter_replace:
    message:
        """
        Filtering alignment {input.sequences} -> {output.sequences}
          - excluding strains in {input.exclude}
          - including strains in {input.include}
          - min length: {params.min_length_query}
        """
    input:
        sequences = "results/{build_name}/masked.fasta",
        metadata = "results/{build_name}/metadata_with_index.tsv",
        # TODO - currently the include / exclude files are not input (origin) specific, but this is possible if we want
        include = "results/{build_name}/include.txt",
        exclude = _collect_exclusion_files,
    output:
        sequences = "results/{build_name}/filtered.fasta",
        filter_log = "results/{build_name}/filtered_log.tsv",
    log:
        "logs/filtered_{build_name}.txt"
    benchmark:
        "benchmarks/filter_{build_name}.txt"
    params:
        min_length_query = _get_filter_min_length_query,
        exclude_where = lambda wildcards: _get_filter_value(wildcards, "exclude_where"),
        min_date = lambda wildcards: _get_filter_value(wildcards, "min_date"),
        ambiguous = lambda wildcards: f"--exclude-ambiguous-dates-by {_get_filter_value(wildcards, 'exclude_ambiguous_dates_by')}" if _get_filter_value(wildcards, "exclude_ambiguous_dates_by") else "",
        date = (date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=500
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            {params.min_length_query} \
            --max-date {params.date} \
            --min-date {params.min_date} \
            {params.ambiguous} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --output {output.sequences} \
            --output-log {output.filter_log} 2>&1 | tee {log};
        """
