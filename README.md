# PolyTrap

**PolyTrap** is a flexible Snakemake pipeline for designing optimized capture probes targeting multiple organisms with customizable priority levels. The pipeline enables researchers to prioritize certain organisms ("focus organisms") while also designing probes for related organisms ("target organisms") with different parameter settings.

## Introduction

PolyTrap was developed to facilitate the design of capture probes for various genomic targets, particularly viruses and other pathogens. By leveraging the power of Snakemake and Conda for workflow management and environment control, PolyTrap offers a reproducible and efficient solution for researchers and diagnosticians aiming to develop robust diagnostic tools. The pipeline is especially valuable for scenarios requiring different design parameters for different target organisms, such as when working with closely related species or strains.

## Features

- **Multi-organism Support:** Design probes for any number of organisms with customizable parameter settings
- **Prioritized Design:** Designate "focus organisms" for enhanced coverage and "target organisms" for standard coverage
- **Enhanced Focus on Segment Flanks:** Optional high-density coverage of segment ends for focus organisms
- **Flexible Configuration:** Organism-specific parameters with defaults for simplicity
- **Blacklist Incorporation:** Exclude unwanted sequences to enhance probe specificity
- **Host Genome Filtering:** Filter candidate probes against host genomes to minimize off-target binding
- **Comprehensive Coverage Analysis:** Evaluate probe coverage across target genomes with organism-specific thresholds
- **Detailed Logging:** Facilitates easy troubleshooting and monitoring of the pipeline
- **Conda Environment Management:** Ensures reproducibility by managing dependencies through Conda

## Requirements

- **Conda:** For environment management. [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- **Snakemake:** Workflow management system.
  ```bash
  conda install -c conda-forge -c bioconda snakemake
  ```

## Installation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yourusername/PolyTrap.git
   cd PolyTrap
   ```

2. **Set Up Conda Environments:**
   
   PolyTrap utilizes Conda environments for managing dependencies. Ensure that the `envs/` directory contains the necessary environment YAML files (`catch.yaml`, `seqkit.yaml`, and `minimap2.yaml`).

   ```bash
   # The environments will be automatically created by Snakemake
   ```

## Configuration

PolyTrap uses a `config.yaml` file to manage parameters and file paths. Below is an example configuration:

```yaml
# Focus organisms (high priority, with special handling)
focus_organisms:
  - name: "puumala"  # Organism identifier
    segments:
      L: "data/puumala_L.fasta"
      M: "data/puumala_M.fasta"
      S: "data/puumala_S.fasta"
  # Add more focus organisms as needed

# Target organisms (regular priority)
target_organisms:
  - name: "orthohanta"  # Organism identifier
    segments:
      L: "data/orthohanta_L.fasta"
      M: "data/orthohanta_M.fasta"
      S: "data/orthohanta_S.fasta"
  # Add more target organisms as needed

# Blacklist FASTA files
blacklist_files:
  - "data/blacklist_1.fasta"
  - "data/blacklist_2.fasta"
  # Add more blacklist files as needed

# Per-organism design parameters
organism_params:
  default:  # Default parameters for all organisms
    probe_length: 120
    mismatches: [2, 3]
    lcf_thres: 100
    coverage: 1.0
    probe_stride: [60, 80]
    design_segment_ends: False
    segment_end_stride: 30
    segment_end_size: 500
    overlap: 50
    use_lsh_minhash: False
    lsh_minhash_threshold: 0.6
  
  puumala:  # Override parameters for specific organism
    mismatches: [2, 3]
    lcf_thres: 100
    design_segment_ends: True
    segment_end_stride: 30
    segment_end_size: 500
    overlap: 50
    analysis_mismatches: 2
    analysis_lcf_thres: 100

# Output directories
output_dirs:
  probes: "results/probes"
  segments: "results/segments"
  blacklist: "results/blacklist"
  analysis: "results/analysis"
  mapping: "results/mapping"

# Other parameters
max_cores: 10
```

### Configuration Keys Explained:

- **focus_organisms**: List of organisms that require special handling (high-density probes, segment end coverage)
- **target_organisms**: List of organisms with standard parameter settings
- **blacklist_files**: FASTA files containing sequences to avoid (e.g., host genomes, closely related species)
- **organism_params**:
  - **default**: Base parameters applied to all organisms (unless overridden)
  - **[organism_name]**: Organism-specific parameter overrides
    - **probe_length**: Length of designed probes (in bp)
    - **mismatches**: Allowed mismatches during probe design
    - **lcf_thres**: Longest common factor threshold
    - **coverage**: Target coverage (0.0-1.0)
    - **probe_stride**: Stride between probe starting positions
    - **design_segment_ends**: Whether to design special high-density probes for segment ends
    - **segment_end_stride**: Stride for segment end probes (smaller value = higher density)
    - **segment_end_size**: Size of segment ends to target (in bp)
    - **overlap**: Overlap between segment middle and segment ends (in bp)
    - **analysis_mismatches**: Mismatches to allow during coverage analysis
    - **analysis_lcf_thres**: LCF threshold during coverage analysis
- **minimap2**: Parameters for minimap2 when mapping probes to host genomes

## Usage

PolyTrap is managed via Snakemake, which orchestrates the workflow based on the `Snakefile` and `config.yaml`.

### Running the Entire Pipeline

Execute the following command from the project root directory:

```bash
snakemake --snakefile baitdesign.smk --configfile config.yaml --cores 10 --use-conda
```

- `--snakefile`: Specifies the path to the Snakefile
- `--configfile`: Specifies the path to the configuration file
- `--cores`: Number of CPU cores to utilize (adjust based on your system)
- `--use-conda`: Enables Conda environment management as specified in the Snakefile

### Dry Run

To perform a dry run and see the jobs that Snakemake would execute without actually running them:

```bash
snakemake --snakefile baitdesign.smk --configfile config.yaml --cores 10 --use-conda -n
```

### Running Specific Rules

If you wish to run only specific rules, such as the final `analyze_probe_coverage` rule, specify its output as the target:

```bash
snakemake --snakefile baitdesign.smk --configfile config.yaml --cores 10 results/analysis/coverage_analysis_2.txt --use-conda
```

This command will execute all necessary rules leading up to the creation of `coverage_analysis_2.txt` for mismatch level 2.

## Workflow Steps

PolyTrap follows these main steps in the workflow:

1. **Extract Segment Ends**: For focus organisms with `design_segment_ends: True`, extract the ends of each segment for high-density probe design
2. **Design Probes**: Design probes for all organisms and segments with varying parameters
3. **Map Probes to Host Genomes**: Map candidate probes to host genomes to identify potential off-target binding
4. **Filter Probes by Host Mapping**: Remove probes that map too frequently to host genomes
5. **Combine and Filter Probes**: Combine all candidate probes and filter to create a non-redundant set
6. **Create Final Probe Sets**: Generate final probe sets for each mismatch level
7. **Analyze Coverage**: Evaluate coverage of the final probe sets against all target genomes

## Troubleshooting

### Common Issues

1. **Conda Environment Errors:**
   - **Solution:** Ensure that Conda is properly installed and that the `envs/` directory contains the necessary YAML files. Use `--use-conda` flag when running Snakemake.

2. **Missing Input Files:**
   - **Solution:** Verify that all input FASTA files are correctly placed in the specified directories and that their paths are correctly specified in `config.yaml`.

3. **Permission Issues:**
   - **Solution:** Ensure you have the necessary read/write permissions for all directories and files involved in the pipeline.

4. **Tool-Specific Errors:**
   - **Solution:** Check the corresponding log files in the `logs/` directory for detailed error messages. Address issues as per the error logs.

### Checking Logs

Each rule writes logs to the `logs/` directory. For example, to view the log for designing probes for a specific organism and segment:

```bash
cat logs/design_probes_regular_puumala_L_ps60_m2.log
```

## Contributing

Contributions are welcome! If you encounter bugs, have feature requests, or wish to contribute code, please follow these steps:

1. **Fork the Repository**
2. **Create a New Branch**
   ```bash
   git checkout -b feature/YourFeatureName
   ```
3. **Commit Your Changes**
   ```bash
   git commit -m "Add your message here"
   ```
4. **Push to the Branch**
   ```bash
   git push origin feature/YourFeatureName
   ```
5. **Open a Pull Request**

Please ensure that your contributions adhere to the project's coding standards and include relevant tests where applicable.

## License

This project is licensed under the [MIT License](LICENSE).
