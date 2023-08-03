ruleorder: filter_replace > filter
ruleorder: annotate_metadata > annotate_metadata_with_index
ruleorder: export_replace > export

rule add_metadata:
    message:
        """
        Deletion to metadata
        """
    input:
        metadata = "nextstrain-ncov-private/metadata.tsv.gz",
        deletions = "orf8ko_profiles/gisaid.washington_ko.tsv"
    params:
        key = "strain",
        add = "ORF8_ko",
        notblank = "coverage"
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
         | tsv-filter \
         -H --not-blank {params.notblank} | gzip > {output.combined}
        """

rule annotate_metadata:
    input:
        metadata="results/{build_name}/metadata_with_nextclade_qc.tsv",
        sequence_index = "results/{build_name}/sequence_index.tsv",
    output:
        metadata="results/{build_name}/metadata_with_index.tsv",
    log:
        "logs/annotate_metadata_with_index_{build_name}.txt",
    benchmark:
        "benchmarks/annotate_metadata_with_index_{build_name}.txt",
    conda: config["conda_environment"]
    shell:
        """
        python3 orf8ko_profiles/annotate_metadata_with_index_add_WA.py \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --output {output.metadata}
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
        tsv-filter -H --str-eq ORF8_ko:Yes --ge coverage:0.95 {input.metadata} | cut -f 1 | sed -n '1!p' > {output.include}
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
        metadata = "results/{build_name}/metadata_filtered.tsv.xz",
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
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --output-log {output.filter_log} 2>&1 | tee {log};
        """

rule find_ko:
    message:
        """
        Calling KO's for SARS-CoV-2 genes
        """
    input:
        align = "results/{build_name}/filtered.fasta"
    params:
        reference = "defaults/reference_seq.gb"
    output:
        ko = "results/{build_name}/ko.tsv"
    shell:
        """
        python orf8ko_profiles/find_ko.py \
            --align {input.align} \
            --ref {params.reference} \
            --output {output.ko}
        """

rule find_ko_with_amplicon:
    message:
        """
        Calling KO's for SARS-CoV-2 genes, including those that fall in known
        amplicon dropout sites
        """
    input:
        align = "results/{build_name}/filtered.fasta"
    params:
        reference = "defaults/reference_seq.gb"
    output:
        ko = "results/{build_name}/ko_with_amplicon.tsv"
    shell:
        """
        python orf8ko_profiles/find_ko.py \
            --align {input.align} \
            --ref {params.reference} \
            --amplicon \
            --output {output.ko}
        """

rule call_clusters:
    message:
        """
        Calling clusters for each gene
        """
    input:
        tree = "results/{build_name}/tree.nwk",
        ko = "results/{build_name}/ko.tsv"
    output:
        dir = directory("results/{build_name}/clusters/"),
        done = touch("results/{build_name}/clusters_done.txt")
    shell:
        """
        python orf8ko_profiles/callClusters.py \
            --tree {input.tree} \
            --ko {input.ko} \
            --outdir {output.dir}
        """

rule call_clusters_with_amplicon:
    message:
        """
        Calling clusters for each gene including KO's that fall in known
        amplicon dropout sites
        """
    input:
        tree = "results/{build_name}/tree.nwk",
        ko = "results/{build_name}/ko_with_amplicon.tsv"
    output:
        dir = directory("results/{build_name}/clusters_with_amplicon/"),
	    done = touch("results/{build_name}/clusters_with_amplicon_done.txt")
    shell:
        """
        python orf8ko_profiles/callClusters.py \
            --tree {input.tree} \
            --ko {input.ko} \
            --outdir {output.dir}
        """

rule add_orf8ko_clusters:
    message:
        """
        Appending orf8 KO clusters
        """
    input:
        metadata = "results/{build_name}/metadata_adjusted.tsv.xz",
        cluster_done = "results/{build_name}/clusters_done.txt",
        #clusters = "results/{build_name}/clusters/clusters_ORF8.tsv",
        ko = "results/{build_name}/ko.tsv"
    params:
        key = "strain",
        add_clusters = "cluster,depth,ORF8_misStops,ORF8_deletions,ORF8_koType",
        add_ko = "ORF8_ko",
        clusters = "results/{build_name}/clusters/clusters_ORF8.tsv"
    output:
        combined = "results/{build_name}/metadata_adjusted_added.tsv.xz"
    shell:
        """
        xz -d -c {input.metadata} | tsv-join \
         --filter-file {params.clusters} \
         --key-fields {params.key} \
         --append-fields {params.add_clusters} \
         -H \
         --write-all ''\
         | tsv-join \
         --filter-file {input.ko} \
          --key-fields {params.key} \
          --append-fields {params.add_ko} \
          --prefix complete \
          -H \
          --write-all ''\
          | xz -c -z > {output.combined}
        """

rule export_replace:
    message: "Exporting data files for Auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/{build_name}/metadata_adjusted_added.tsv.xz",
        node_data = _get_node_data_by_wildcards,
        auspice_config = get_auspice_config,
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"].get(w.build_name, {}) else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"],
        description = rules.build_description.output.description
    output:
        auspice_json = "results/{build_name}/ncov_with_accessions.json",
        root_sequence_json = "results/{build_name}/ncov_with_accessions_root-sequence.json"
    log:
        "logs/export_{build_name}.txt"
    benchmark:
        "benchmarks/export_{build_name}.txt"
    params:
        title = export_title
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --minify-json \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
