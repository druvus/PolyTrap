import os
import glob
from itertools import product

configfile: "config.yaml"

################################################
# 1) Load Configuration
################################################

# Load focus and target organisms
focus_organisms = config.get("focus_organisms", [])
target_organisms = config.get("target_organisms", [])

# Load blacklist files
blacklist_files = config.get("blacklist_files", [])

# Load default organism parameters
default_params = config.get("organism_params", {}).get("default", {})

# Output directories and other settings
output_dirs = config["output_dirs"]
conda_envs = config["conda_envs"]
logs_dir = config["logs_dir"]

# Create logs directory
os.makedirs(logs_dir, exist_ok=True)

################################################
# 2) Utility Functions
################################################

def resolve_organism_params(organism_name):
    """
    Resolves parameters for an organism by merging default params with organism-specific overrides.
    """
    # Start with default parameters
    params = dict(default_params)
    
    # Override with organism-specific parameters if available
    organism_specific = config.get("organism_params", {}).get(organism_name, {})
    params.update(organism_specific)
    
    return params

def get_all_organisms():
    """
    Returns a list of all organism names.
    """
    focus_names = [org["name"] for org in focus_organisms]
    target_names = [org["name"] for org in target_organisms]
    return focus_names + target_names

def get_all_segments():
    """
    Returns a dict mapping organism names to their segments.
    """
    segments = {}
    
    # Add focus organism segments
    for org in focus_organisms:
        segments[org["name"]] = list(org["segments"].keys())
    
    # Add target organism segments
    for org in target_organisms:
        segments[org["name"]] = list(org["segments"].keys())
    
    return segments

def get_segment_path(organism_name, segment_name):
    """
    Returns the path to a segment file for a given organism.
    """
    # Check focus organisms
    for org in focus_organisms:
        if org["name"] == organism_name:
            return org["segments"].get(segment_name)
    
    # Check target organisms
    for org in target_organisms:
        if org["name"] == organism_name:
            return org["segments"].get(segment_name)
    
    return None

def is_focus_organism(organism_name):
    """
    Returns True if the organism is a focus organism, False otherwise.
    """
    return organism_name in [org["name"] for org in focus_organisms]

def design_command_params(organism, segment, params, output, input_file):
    """
    Constructs the command for design.py based on the given parameters.
    """
    cmd = [
        "design.py",
        "-pl", str(params.get("probe_length", 120)),
        "-ps", str(params.get("probe_stride", 60)),
        "-m", str(params.get("mismatches", 2)),
        "-l", str(params.get("lcf_thres", 100)),
        "--skip-set-cover"
    ]
    
    for blk in blacklist_files:
        cmd.extend(["--avoid-genomes", blk])
    
    if params.get("use_lsh_minhash", False):
        cmd.extend(["--filter-with-lsh-minhash", str(params.get("lsh_minhash_threshold", 0.6))])
    
    cmd.extend(["-o", output, input_file])
    
    return " ".join(cmd)

def get_all_mismatches():
    """Get unique mismatch values from all organism configs."""
    all_mismatches = set()
    
    for org_name in get_all_organisms():
        params = resolve_organism_params(org_name)
        mismatches = params.get("mismatches", [2])
        
        if isinstance(mismatches, list):
            all_mismatches.update(mismatches)
        else:
            all_mismatches.add(mismatches)
    
    return sorted(list(all_mismatches))

################################################
# 3) Expansions for Designing Probes
################################################

# Get all organisms and their segments
all_organisms = get_all_organisms()
all_segments = get_all_segments()

# Create expansions for all probes
all_probes = []
segment_ends_needed = set()

# For each organism and segment, add regular probes
for org_name in all_organisms:
    params = resolve_organism_params(org_name)
    
    # Add to mismatch and stride combinations
    mismatches = params.get("mismatches", [2])
    if not isinstance(mismatches, list):
        mismatches = [mismatches]
    
    probe_strides = params.get("probe_stride", [60])
    if not isinstance(probe_strides, list):
        probe_strides = [probe_strides]
    
    # For each segment and parameter combination
    for segment in all_segments.get(org_name, []):
        for m, ps in product(mismatches, probe_strides):
            probe_file = os.path.join(
                output_dirs["probes"], 
                f"{org_name}_{segment}_ps{ps}_m{m}_probes.fasta"
            )
            all_probes.append(probe_file)
    
    # If segment ends should be designed separately
    if params.get("design_segment_ends", False):
        end_stride = params.get("segment_end_stride", 30)
        for segment in all_segments.get(org_name, []):
            for m in mismatches:
                probe_file = os.path.join(
                    output_dirs["probes"], 
                    f"{org_name}_{segment}_ends_ps{end_stride}_m{m}_probes.fasta"
                )
                all_probes.append(probe_file)
                # Track which segment ends we need to extract
                segment_ends_needed.add((org_name, segment))

# Expand the segment ends that need to be extracted
segment_ends_files = [
    os.path.join(output_dirs["segments"], f"{org}_{segment}_ends.fasta") 
    for org, segment in segment_ends_needed
]

################################################
# 4) Identify Genome References for Mapping
################################################

genome_folder = config["genome_folder"]
genome_files = glob.glob(os.path.join(genome_folder, "*", "assembly", "*.fa.gz"))

# Build a list of dicts for host, assembly, and file path
genomes = [
    {
        "host": os.path.basename(os.path.dirname(os.path.dirname(filepath))),
        "assembly": os.path.splitext(os.path.basename(filepath))[0],
        "filepath": filepath
    }
    for filepath in genome_files
]

# Unique IDs like "Aedes_aegypti_GCA_002204515.1"
genome_ids = [
    f"{g['host']}_{g['assembly'].replace('.', '_')}" 
    for g in genomes
]

################################################
# 5) Construct expansions for mapping each design
################################################

def design_file_stem(path):
    """
    Given an absolute path like:
        results/probes/puumala_L_ps60_m2_probes.fasta
    Return:
        puumala_L_ps60_m2_probes
    """
    base = os.path.basename(path)       # puumala_L_ps60_m2_probes.fasta
    stem, _ = os.path.splitext(base)    # puumala_L_ps60_m2_probes
    return stem

# For each design FASTA, map to each genome -> a PAF
paf_files = [
    os.path.join(
        output_dirs["mapping"],
        f"{design_file_stem(design_file)}__{genome_id}.paf"  # Note the double underscore
    )
    for design_file in all_probes
    for genome_id in genome_ids
]

# Filtered matching oligos for each PAF
matching_oligos_txt = [
    paf.replace(".paf", "_matching_oligos.txt")
    for paf in paf_files
]

################################################
# 6) rule all
################################################

wildcard_constraints:
    organism="|".join(get_all_organisms()),
    segment="|".join([seg for org_segs in all_segments.values() for seg in org_segs]),
    ps="\d+",
    m="\d+"

rule all:
    """
    The 'all' rule ensures all components of the pipeline are executed.
    """
    input:
        # 1) All design FASTAs
        all_probes,
        # 2 & 3) Mapping outputs (PAF + matching oligos per design)
        paf_files,
        matching_oligos_txt,
        # 4) Combined candidate probes
        os.path.join(output_dirs["probes"], "all_candidate_probes.fasta"),
        # 5) Final sets for each mismatch value
        expand(
            os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta"),
            m=get_all_mismatches()
        ),
        # 6) Coverage analysis for each mismatch value
        expand(
            os.path.join(output_dirs["analysis"], "coverage_analysis_{m}.txt"),
            m=get_all_mismatches()
        ),
        # 7) Host genome mapping analysis
        os.path.join(output_dirs["analysis"], "host_probe_counts.txt"),
        # 8) Filtered host mapping counts for final probes
        os.path.join(output_dirs["analysis"], "final_probes_host_counts.txt")

################################################
# Step 1: Extract segments (for organisms that need segment ends)
################################################

rule extract_segment_ends:
    """
    Extracts the ends and middles from segments for organisms that need segment-end probe design.
    """
    input:
        genome=lambda w: get_segment_path(w.organism, w.segment)
    output:
        ends=os.path.join(output_dirs["segments"], "{organism}_{segment}_ends.fasta"),
        middles=os.path.join(output_dirs["segments"], "{organism}_{segment}_middles.fasta")
    params:
        params=lambda w: resolve_organism_params(w.organism),
        end_size=lambda w: resolve_organism_params(w.organism).get("segment_end_size", 500),
        overlap=lambda w: resolve_organism_params(w.organism).get("overlap", 50),
        start_middle=lambda w: resolve_organism_params(w.organism).get("segment_end_size", 500) + 1 - resolve_organism_params(w.organism).get("overlap", 50),
        end_middle=lambda w: -(resolve_organism_params(w.organism).get("segment_end_size", 500) + 1 - resolve_organism_params(w.organism).get("overlap", 50)),
        output_dir=output_dirs["segments"],
        temp_dir=lambda w: os.path.join(output_dirs["segments"], f"temp_{w.organism}_{w.segment}")
    log:
        os.path.join(logs_dir, "extract_segment_ends_{organism}_{segment}.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        mkdir -p {params.temp_dir}
        seqkit subseq --region 1:{params.end_size} --update-faidx {input.genome} \
          > {params.temp_dir}/{wildcards.organism}_{wildcards.segment}_5prime.fasta 2>> {log}
        seqkit subseq --region -{params.end_size}:-1 --update-faidx {input.genome} \
          > {params.temp_dir}/{wildcards.organism}_{wildcards.segment}_3prime.fasta 2>> {log}

        cat {params.temp_dir}/{wildcards.organism}_{wildcards.segment}_5prime.fasta \
            {params.temp_dir}/{wildcards.organism}_{wildcards.segment}_3prime.fasta \
          > {output.ends} 2>> {log}

        seqkit subseq --region {params.start_middle}:{params.end_middle} --update-faidx {input.genome} \
          > {output.middles} 2>> {log}

        rm -r {params.temp_dir}
        """

################################################
# Step 2: Design Probes
################################################

rule design_probes_segment_ends:
    """
    Designs probes specifically for segment ends, if enabled for the organism.
    """
    input:
        ends=os.path.join(output_dirs["segments"], "{organism}_{segment}_ends.fasta")
    output:
        probes=os.path.join(output_dirs["probes"], "{organism}_{segment}_ends_ps{ps}_m{m}_probes.fasta")
    params:
        organism_params=lambda w: resolve_organism_params(w.organism),
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: (
            f"touch {output.probes}" if not resolve_organism_params(wildcards.organism).get("design_segment_ends", False)
            else design_command_params(
                wildcards.organism,
                wildcards.segment,
                {
                    "probe_length": resolve_organism_params(wildcards.organism).get("probe_length", 120),
                    "probe_stride": int(wildcards.ps),
                    "mismatches": int(wildcards.m),
                    "lcf_thres": resolve_organism_params(wildcards.organism).get("lcf_thres", 100),
                },
                output.probes,
                input.ends
            )
        )
    log:
        os.path.join(logs_dir, "design_probes_segment_ends_{organism}_{segment}_ps{ps}_m{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

rule design_probes_regular:
    """
    Designs probes for all organism segments with various stride/mismatch combos.
    """
    input:
        genome=lambda w: get_segment_path(w.organism, w.segment)
    output:
        probes=os.path.join(output_dirs["probes"], "{organism}_{segment}_ps{ps}_m{m}_probes.fasta")
    params:
        organism_params=lambda w: resolve_organism_params(w.organism),
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: design_command_params(
            wildcards.organism,
            wildcards.segment,
            {
                "probe_length": resolve_organism_params(wildcards.organism).get("probe_length", 120),
                "probe_stride": int(wildcards.ps),
                "mismatches": int(wildcards.m),
                "lcf_thres": resolve_organism_params(wildcards.organism).get("lcf_thres", 100),
            },
            output.probes,
            input.genome
        )
    log:
        os.path.join(logs_dir, "design_probes_regular_{organism}_{segment}_ps{ps}_m{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

################################################
# Step 3: Map Each Design FASTA to Each Genome
################################################
rule map_probes_to_genome:
    """
    Maps each designed probe FASTA to each genome reference.
    """
    input:
        probes=lambda wc: next(
            path for path in all_probes
            if design_file_stem(path) == wc.design_file_stem
        ),
        genome=lambda wc: next(
            g["filepath"] for g in genomes
            if f"{g['host']}_{g['assembly'].replace('.', '_')}" == wc.genome_id
        )
    output:
        paf=os.path.join(output_dirs["mapping"], "{design_file_stem}__{genome_id}.paf")
    params:
        preset=config["minimap2"]["preset"],
        k=config["minimap2"]["k"],
        w=config["minimap2"]["w"],
        threads=config["minimap2"]["threads"],
        output_dir=output_dirs["mapping"]
    log:
        os.path.join(logs_dir, "map_probes_to_genome_{design_file_stem}__{genome_id}.log")
    conda:
        os.path.join(conda_envs, "minimap2.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        minimap2 -x {params.preset} -k {params.k} -w {params.w} -t {params.threads} \
          {input.genome} {input.probes} \
          > {output.paf} 2> {log}
        """

rule filter_paf:
    """
    Applies awk to each PAF file to select records meeting coverage/identity thresholds.
    """
    input:
        paf=os.path.join(output_dirs["mapping"], "{design_file_stem}__{genome_id}.paf")
    output:
        matching_oligos=os.path.join(
            output_dirs["mapping"],
            "{design_file_stem}__{genome_id}_matching_oligos.txt"
        )
    params:
        awk_condition=f"'$1 >= {config['minimap2']['min_coverage']} && $10/$11 >= {config['minimap2']['min_identity']}'"
    log:
        os.path.join(logs_dir, "filter_paf_{design_file_stem}__{genome_id}.log")
    shell:
        """
        awk {params.awk_condition} {input.paf} | cut -f1 | sort -u \
          > {output.matching_oligos} 2> {log}
        """

rule count_host_mappings:
    """
    Counts how many times each oligo is mapped across all host genomes.
    """
    input:
        matching_oligos=matching_oligos_txt
    output:
        counts=os.path.join(output_dirs["analysis"], "host_probe_counts.txt")
    params:
        output_dir=output_dirs["analysis"]
    log:
        os.path.join(logs_dir, "count_host_mappings.log")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.matching_oligos} | sort | uniq -c | sort -nr > {output.counts} 2> {log}
        """

rule filter_high_host_hits:
    input:
        probes=os.path.join(output_dirs["probes"], "{organism}_{segment}_{suffix}_probes.fasta"),
        host_counts=os.path.join(output_dirs["analysis"], "host_probe_counts.txt") 
    output:
        filtered=os.path.join(output_dirs["probes"], "{organism}_{segment}_{suffix}_host_filtered_probes.fasta")
    params:
        max_hits=lambda w: config.get("analysis", {}).get("max_host_hits", 2),
        temp_exclude=lambda w: f"temp_exclude_{w.organism}_{w.segment}_{w.suffix}.txt"
    wildcard_constraints:
        suffix="(ps\d+_m\d+|ends_ps\d+_m\d+)"
    log:
        os.path.join(logs_dir, "filter_high_host_hits_{organism}_{segment}_{suffix}.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        if [ -s {input.probes} ]; then
            # Create exclude list with unique temp file
            awk '$1 > {params.max_hits} {{print substr($0,9)}}' {input.host_counts} > {params.temp_exclude}
            
            if [ -s {params.temp_exclude} ]; then
                seqkit grep -v -f {params.temp_exclude} {input.probes} -w 0 > {output.filtered} 2> {log}
            else
                cp {input.probes} {output.filtered}
            fi
            rm -f {params.temp_exclude}
        else
            touch {output.filtered}
        fi
        """

################################################
# Step 4: Combine + Filter
################################################

rule combine_candidate_probes:
    input:
        probes=expand(
            os.path.join(output_dirs["probes"], "{design}_host_filtered_probes.fasta"),
            design=[p.replace("_probes.fasta", "") for p in map(os.path.basename, all_probes)]
        )
    output:
        combined=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta")
    log:
        os.path.join(logs_dir, "combine_candidate_probes.log")
    shell:
        "cat {input.probes} > {output.combined}"

rule filter_probes:
    """
    Filters the combined candidate probes using references from all organisms.
    """
    input:
        candidate_probes=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta"),
        segment_genomes=lambda w: [
            get_segment_path(org, seg) 
            for org in get_all_organisms() 
            for seg in all_segments.get(org, [])
        ],
        segment_ends=segment_ends_files
    output:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes_{m}.fasta")
    params:
        pl=default_params.get("probe_length", 120),
        l=min([resolve_organism_params(org).get("lcf_thres", 100) for org in get_all_organisms()]),
        c=1.0,
        output_dir=output_dirs["probes"],
        all_references=lambda w, input: " ".join(input.segment_genomes + input.segment_ends)
    wildcard_constraints:
        m="|".join(map(str, get_all_mismatches()))
    log:
        os.path.join(logs_dir, "filter_probes_{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        design.py \
            {params.all_references} \
            -pl {params.pl} \
            -m {wildcards.m} \
            -l {params.l} \
            -c {params.c} \
            --filter-from-fasta {input.candidate_probes} \
            -o {output.filtered_probes} > {log} 2>&1
        """

rule final_probe_set:
    """
    Creates final probe sets for each mismatch value.
    """
    input:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes_{m}.fasta")
    output:
        final=os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta")
    params:
        output_dir=output_dirs["probes"]
    wildcard_constraints:
        m="|".join(map(str, get_all_mismatches()))
    shell:
        """
        mkdir -p {params.output_dir}
        cp {input.filtered_probes} {output.final}
        """

rule analyze_host_matches_in_final:
    """
    Identifies which oligos from host_probe_counts appear in final probe sets.
    """
    input:
        host_counts=os.path.join(output_dirs["analysis"], "host_probe_counts.txt"),
        final_probes=expand(
            os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta"),
            m=get_all_mismatches()
        )
    output:
        filtered_counts=os.path.join(output_dirs["analysis"], "final_probes_host_counts.txt")
    log:
        os.path.join(logs_dir, "analyze_host_matches_in_final.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        # Extract sequence IDs from final probe sets
        cat {input.final_probes} | grep "^>" | cut -d ">" -f 2 > temp_probe_ids.txt

        # Extract oligos and their counts, keeping only those in final probe sets
        awk 'NR==FNR {{probes[$1]=1; next}} 
             {{oligo=substr($0,9); if(probes[oligo]) print $0}}' \
             temp_probe_ids.txt {input.host_counts} > {output.filtered_counts} 2> {log}
        
        rm temp_probe_ids.txt
        """

################################################
# Step 5: Validate (Coverage Analysis)
################################################

rule analyze_probe_coverage:
    """
    Uses analyze_probe_coverage.py to perform coverage analysis for each mismatch value.
    """
    input:
        probes=os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta"),
        all_segments=lambda w: [
            get_segment_path(org, seg)
            for org in get_all_organisms()
            for seg in all_segments.get(org, [])
        ],
        segment_ends=segment_ends_files
    output:
        coverage=os.path.join(output_dirs["analysis"], "coverage_analysis_{m}.txt")
    params:
        output_dir=output_dirs["analysis"],
        focus_segments=lambda w: " ".join([
            get_segment_path(org, seg)
            for org in [o["name"] for o in focus_organisms]
            for seg in all_segments.get(org, [])
        ]),
        target_segments=lambda w: " ".join([
            get_segment_path(org, seg)
            for org in [o["name"] for o in target_organisms]
            for seg in all_segments.get(org, [])
        ]),
        focus_ends=lambda w: " ".join([
            os.path.join(output_dirs["segments"], f"{org}_{seg}_ends.fasta")
            for org in [o["name"] for o in focus_organisms]
            for seg in all_segments.get(org, [])
            if (org, seg) in segment_ends_needed
        ]),
        focus_params=lambda w: {
            "m": w.m,
            "l": min([
                resolve_organism_params(org["name"]).get("analysis_lcf_thres",
                resolve_organism_params(org["name"]).get("lcf_thres", 100))
                for org in focus_organisms
            ])
        },
        target_params=lambda w: {
            "m": w.m, 
            "l": min([
                resolve_organism_params(org["name"]).get("analysis_lcf_thres",
                resolve_organism_params(org["name"]).get("lcf_thres", 90))
                for org in target_organisms
            ])
        }
    wildcard_constraints:
        m="|".join(map(str, get_all_mismatches()))
    log:
        os.path.join(logs_dir, "analyze_probe_coverage_{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        echo "=== Coverage Analysis Report - $(date) ===" >> {log}
        
        echo "Analyzing coverage for focus organisms full genomes..." >> {log}
        analyze_probe_coverage.py \
            -d {params.focus_segments} \
            -f {input.probes} \
            -m {params.focus_params[m]} \
            -l {params.focus_params[l]} \
            --print-analysis \
          > {params.output_dir}/focus_full_coverage_{wildcards.m}.txt 2>> {log}
        echo "Completed coverage analysis for focus organisms full genomes." >> {log}
        
        if [ ! -z "{params.focus_ends}" ]; then
          echo "\nAnalyzing coverage for focus organisms flanks..." >> {log}
          analyze_probe_coverage.py \
              -d {params.focus_ends} \
              -f {input.probes} \
              -m {params.focus_params[m]} \
              -l {params.focus_params[l]} \
              --print-analysis \
            > {params.output_dir}/focus_flanks_coverage_{wildcards.m}.txt 2>> {log}
          echo "Completed coverage analysis for focus organisms flanks." >> {log}
        else
          touch {params.output_dir}/focus_flanks_coverage_{wildcards.m}.txt
        fi
        
        echo "\nAnalyzing coverage for target organisms..." >> {log}
        analyze_probe_coverage.py \
            -d {params.target_segments} \
            -f {input.probes} \
            -m {params.target_params[m]} \
            -l {params.target_params[l]} \
            --print-analysis \
          > {params.output_dir}/target_coverage_{wildcards.m}.txt 2>> {log}
        echo "Completed coverage analysis for target organisms." >> {log}
        
        echo "\nCombining all coverage reports into a single file..." >> {log}
        cat {params.output_dir}/focus_full_coverage_{wildcards.m}.txt \
            {params.output_dir}/focus_flanks_coverage_{wildcards.m}.txt \
            {params.output_dir}/target_coverage_{wildcards.m}.txt \
          > {output.coverage} 2>> {log}
        echo "Combined coverage analysis report created at {output.coverage}" >> {log}
        echo "=== Coverage Analysis Completed ===" >> {log}
        """
