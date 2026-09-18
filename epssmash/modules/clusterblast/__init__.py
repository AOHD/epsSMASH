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

from antismash.common import path, json
from antismash.common.secmet import Record, Region
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs
from antismash.modules.clusterblast import (
    ClusterBlastResults,
    check_clusterblast_files,
    check_options,
    regenerate_previous_results,
    get_result_limit,
    load_clusterblast_database,
    prepare_single_database,
    regenerate_previous_results,
    will_handle,
)
from antismash.modules.clusterblast.core import (
    get_core_gene_ids,
    load_reference_clusters_from_dir,
    load_reference_proteins_from_dir,
    parse_all_clusters,
    run_diamond_on_all_regions,
    score_clusterblast_output,
)
from antismash.modules.clusterblast.data_structures import (
    Protein,
    ReferenceCluster,
)
from antismash.modules.clusterblast.results import GeneralResults, RegionResult
from antismash.modules.clusterblast.html_output import (
    generate_div,
    generate_javascript_data as _original_gen_js_data,
)

NAME = "clusterblast"
SHORT_DESCRIPTION = "Runs clusterblast over custom data"


def generate_html(region_layer: RegionLayer, results: ClusterBlastResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    html = HTMLSections("clusterblast")
    region = region_layer.region_feature

    references = [ref for ref, _ in results.general.region_results[region.get_region_number() - 1].ranking]

    base_tooltip = ("Shows %s that are similar to the current region. Genes marked with the "
                    "same colour are interrelated. White genes have no relationship.<br>"
                    "Click on reference genes to show details of similarities to "
                    "genes within the current region.")

    if options_layer.cb_general or region.clusterblast is not None:
        tooltip = base_tooltip % "regions from the epsSMASH database of manually validated EPS gene clusters"
        #tooltip += "<br>Click on an accession to open that entry in the antiSMASH database (if applicable)."
        div = generate_div(region_layer, record_layer, options_layer, "clusterblast",
                           tooltip, references, title="Similar gene clusters")
        html.add_detail_section("Clusterblast", div, "clusterblast")

    return html


def generate_javascript_data(record: Record, region: Region, results: ClusterBlastResults,
                             ) -> json.JSONBase:
    data = _original_gen_js_data(record, region, results)
    # reusing the name "clusterblast" causes the original result generation to
    # including linking to the antiSMASH-DB, so disable URL(s) in the output
    for variant in data.references:
        variant.url = ""
    return data


def get_arguments() -> ModuleArgs:
    """ Builds the args for the clusterblast module """
    args = ModuleArgs('Clusterblast options', 'cb')
    args.add_analysis_toggle('general',
                             dest='general',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against a "
                                  "database of epsSMASH-predicted clusters.")
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
    return options.cb_general


def load_reference_clusters(searchtype: str) -> dict[str, ReferenceCluster]:  # pylint: disable=unused-argument
    """ Load gene cluster database

        Arguments:
            searchtype: determines which database to use, allowable values

        Returns:
            a dictionary mapping reference cluster name to ReferenceCluster
            instance
    """
    options = get_config()

    logging.info("Clusterblast: Loading gene cluster database into memory...")
    if options.database_dir is None:
        raise ValueError("No database directory specified")
    data_dir = os.path.join(options.database_dir, 'clusterblast')

    return load_reference_clusters_from_dir(data_dir)


def load_reference_proteins(searchtype: str) -> dict[str, Protein]:
    """ Load protein database

        Arguments:
            searchtype: determines which database to use

        Returns:
            a dictionary mapping protein name to Protein instance
    """
    options = get_config()

    logging.info("ClusterBlast: Loading gene cluster database proteins into memory...")
    data_dir = os.path.join(options.database_dir, 'clusterblast')

    return load_reference_proteins_from_dir(data_dir)


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
    database = os.path.join(options.database_dir, 'clusterblast', 'proteins.fasta')
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
    failure_messages: list[str] = []

    # general
    clusterblastdir = os.path.join(get_config().database_dir, "clusterblast")
    failure_messages.extend(prepare_single_database(clusterblastdir,
                                                    raise_exceptions=not logging_only))
    return failure_messages


def run_on_record(record: Record, results: Optional[ClusterBlastResults],
                  options: ConfigType) -> ClusterBlastResults:
    """ Runs over the given record and finds similar areas in the databases """
    if not results:
        results = ClusterBlastResults(record.id)
    if options.cb_general and not results.general:
        clusters, proteins = load_clusterblast_database()
        results.general = perform_clusterblast(options, record, clusters, proteins)
    return results
