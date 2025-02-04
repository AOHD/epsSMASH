# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
An override of antiSMASH's clusterblast module, excluding unnecessarsy/unused portions.
"""

import logging
import os
from typing import Optional

import antismash
from antismash.common.html_renderer import HTMLSections
from antismash.common.layers import (
    RecordLayer,
    RegionLayer,
    OptionsLayer,
)

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs
from antismash.modules.clusterblast import (
    ClusterBlastResults,
    check_clusterblast_files,
    check_options,
    get_result_limit,
    load_clusterblast_database,
    will_handle,
    run_knownclusterblast_on_record,
    prepare_known_data,
)
from antismash.modules.clusterblast.core import (
    get_core_gene_ids,
    parse_all_clusters,
    run_diamond_on_all_regions,
    score_clusterblast_output,
)
from antismash.modules.clusterblast.data_structures import (
    Protein,
    ReferenceCluster,
)
from antismash.modules.clusterblast.results import GeneralResults, RegionResult
from antismash.modules.clusterblast.html_output import generate_div

NAME = "customblast"
SHORT_DESCRIPTION = "Runs clusterblast over custom data"


def generate_html(region_layer: RegionLayer, results: ClusterBlastResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    html = HTMLSections("clusterblast")
    region = region_layer.region_feature

    base_tooltip = ("Shows %s that are similar to the current region. Genes marked with the "
                    "same colour are interrelated. White genes have no relationship.<br>"
                    "Click on reference genes to show details of similarities to "
                    "genes within the current region.")

    if options_layer.cb_general or region.clusterblast is not None:
        tooltip = base_tooltip % "regions from the antiSMASH database"
        tooltip += "<br>Click on an accession to open that entry in the antiSMASH database (if applicable)."
        div = generate_div(region_layer, record_layer, options_layer, "clusterblast", tooltip)
        html.add_detail_section("Customblast", div, "clusterblast")

    return html


def get_arguments() -> ModuleArgs:
    """ Builds the args for the clusterblast module """
    args = ModuleArgs('CustomBlast options', 'cb')
    args.add_analysis_toggle('general',
                             dest='general',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against a "
                                  "database of epsSMASH-predicted clusters.")
    args.add_analysis_toggle('knownclusters',
                             dest='knownclusters',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against known "
                                  "gene clusters from the literature.")
    args.add_option('min-homology-scale',
                    dest='min_homology_scale',
                    metavar="LIMIT",
                    type=float,
                    default=0.0,
                    help="A minimum scaling factor for the query BGC in ClusterBlast results."
                         " Valid range: 0.0 - 1.0.   "
                         " Warning: some homologous genes may no longer be visible!"
                         " (default: %(default)s)")
    args.add_option('nclusters',
                    dest='nclusters',
                    metavar="count",
                    type=int,
                    default=10,
                    help="Number of clusters to display,"
                         f" cannot be greater than {get_result_limit()}. (default: %(default)s)")
    return args


def check_prereqs(options: ConfigType) -> list[str]:
    "Check if all required applications are around"
    _required_binaries = [
        'blastp',
        'makeblastdb',
        'diamond'
    ]

    failure_messages = []
    for binary_name in _required_binaries:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    if "diamond" not in get_config().executables:
        failure_messages.append("cannot check clusterblast databases, no diamond executable present")
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    return options.cb_general or options.cb_knownclusters


def load_reference_clusters(searchtype: str) -> dict[str, ReferenceCluster]:
    """ Load gene cluster database

        Arguments:
            searchtype: determines which database to use, allowable values:
                            clusterblast

        Returns:
            a dictionary mapping reference cluster name to ReferenceCluster
            instance
    """
    options = get_config()

    if searchtype == "clusterblast":
        logging.info("CustomBlast: Loading gene cluster database into memory...")
        if options.database_dir is None:
            raise ValueError("No database directory specified")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database into memory...")
        kcb_root = os.path.join(options.database_dir, "knownclusterblast")
        version = path.find_latest_database_version(kcb_root)
        data_dir = os.path.join(kcb_root, version)
        print(data_dir)    

    reference_cluster_file = os.path.join(data_dir, "clusters.txt")
    with open(reference_cluster_file, "r", encoding="utf-8") as handle:
        filetext = handle.read()
    lines = [line for line in filetext.splitlines() if "\t" in line]
    clusters = {}
    for i in lines:
        tabs = i.split("\t")
        accession = tabs[0]
        description = tabs[1]
        cluster_number = tabs[2]
        cluster_type = tabs[3]
        tags = tabs[4].split(";")
        proteins = tabs[5].split(";")
        if not proteins[-1]:
            proteins.pop(-1)
        cluster = ReferenceCluster(accession, cluster_number, proteins,
                                   description, cluster_type, tags)
        clusters[cluster.get_name()] = cluster
    return clusters


def load_reference_proteins(searchtype: str) -> dict[str, Protein]:
    """ Load protein database

        Arguments:
            searchtype: determines which database to use, allowable values:
                            clusterblast, subclusterblast, knownclusterblast
        Returns:
            a dictionary mapping protein name to Protein instance
    """
    options = get_config()
    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = os.path.join(options.database_dir, 'customblast')
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database proteins into memory...")
        kcb_root = os.path.join(options.database_dir, "knownclusterblast")
        version = path.find_latest_database_version(kcb_root)
        data_dir = os.path.join(kcb_root, version)
        
    
    protein_file = os.path.join(data_dir, "proteins.fasta")
    proteins = {}
    with open(protein_file, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line[0] != ">":
                continue
            # some lines are malformed, so always split the name off the annotation
            # e.g. >x|y|1-2|-|z|Urea_carboxylase_{ECO:0000313|EMBL:CCF11062.1}|CRH36422
            tabs = line.split("|", 5)
            annotations, name = tabs[5].rsplit("|", 1)
            locustag = tabs[4]
            location = tabs[2]
            strand = tabs[3]
            proteins[locustag] = Protein(name, locustag, location, strand, annotations)
    return proteins


antismash.modules.clusterblast.core.load_reference_clusters = load_reference_clusters
antismash.modules.clusterblast.core.load_reference_proteins = load_reference_proteins


def perform_clusterblast(options: ConfigType, record: Record,
                         db_clusters: dict[str, ReferenceCluster],
                         db_proteins: dict[str, Protein]) -> GeneralResults:
    """ Run BLAST on gene cluster proteins for each cluster, parse output and
        return result rankings for each cluster

        Arguments:
            options: antismash Config
            record: the Record to analyse
            db_clusters: a dict mapping reference cluster name to ReferenceCluster
            db_proteins: a dict mapping reference protein name to Protein

        Returns:
            a GeneralResults instance with results for each cluster in the record
    """
    regions = record.get_regions()
    database = os.path.join(options.database_dir, 'customblast', 'proteins.fasta')
    blastoutput = run_diamond_on_all_regions(regions, database)

    clusters_by_number, _ = parse_all_clusters(blastoutput, record,
                                               min_seq_coverage=10,
                                               min_perc_identity=30)
    results = GeneralResults(record.id)

    core_gene_accessions = get_core_gene_ids(record)

    for region in regions:
        region_number = region.get_region_number()
        cluster_names_to_queries = clusters_by_number.get(region_number, {})
        ranking = score_clusterblast_output(db_clusters, core_gene_accessions,
                                            cluster_names_to_queries)
        # store the results
        result = RegionResult(region, ranking, db_proteins, "general")

        results.add_region_result(result, db_clusters, db_proteins)

    results.write_to_file(record, options)
    return results


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Prepare the databases. """
    failure_messages = []
    # known
    failure_messages.extend(prepare_known_data(logging_only))

    # general
    clusterblastdir = os.path.join(get_config().database_dir, "clusterblast")
    if "mounted_at_runtime" in clusterblastdir:  # can't prepare these
        return failure_messages
    cluster_defs = os.path.join(clusterblastdir, 'clusters.txt')
    protein_seqs = os.path.join(clusterblastdir, "proteins.fasta")
    db_file = os.path.join(clusterblastdir, "proteins.dmnd")

    # check the DBv3 region info exists instead of single cluster numbers
    with open(protein_seqs, encoding="utf-8") as handle:
        sample = handle.readline()
    if "-" not in sample.split("|", 3)[1]:
        failure_messages.append("clusterblast database out of date, update with download-databases")
        # and don't bother pressing them
        return failure_messages

    failure_messages.extend(check_clusterblast_files(cluster_defs, protein_seqs, db_file, logging_only=logging_only))

    return failure_messages


def run_on_record(record: Record, results: Optional[ClusterBlastResults],
                  options: ConfigType) -> ClusterBlastResults:
    """ Runs over the given record and finds similar areas in the databases """
    if not results:
        results = ClusterBlastResults(record.id)
    if options.cb_general and not results.general:
        clusters, proteins = load_clusterblast_database()
        results.general = perform_clusterblast(options, record, clusters, proteins)
    if options.cb_knownclusters and not results.knowncluster:
        results.knowncluster = run_knownclusterblast_on_record(record, options)
    return results
