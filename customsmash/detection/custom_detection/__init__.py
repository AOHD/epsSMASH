# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
A direct replacement of antiSMASH's cluster detection module, using all the same logic,
but with custom rules and profiles. There's no need to reuse any or all of that module,
as long as the general rules for inputs and outputs of antiSMASH modules are respected.

Some functions are expected by antiSMASH, for details antismash.custom_typing.pyi

This module could instead be implemented as a class, using the above information and
inheriting antismash.custom_typing.AntismashModule
"""
import logging
from typing import Any, Optional, Self

# import any components being reused from antiSMASH
from antismash.common import path
from antismash.common.hmmer import ensure_database_pressed
from antismash.common.hmm_rule_parser.cluster_prediction import (
    create_rules,
    detect_protoclusters_and_signatures,
    RuleDetectionResults,
    Ruleset,
)
from antismash.config import ConfigType
from antismash.common.module_results import DetectionResults
from antismash.common.secmet.record import Record
from antismash.common.secmet.features import Protocluster
from antismash.common.signature import HmmSignature
from antismash.common.signature import get_signature_profiles
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage
from antismash.detection.hmm_detection import check_prereqs as original_check_prereqs

NAME = "customsmash_detection"
SHORT_DESCRIPTION = "some kind of protocluster detection"
# the detection stage defines when the module is run in the detection process
DETECTION_STAGE = DetectionStage.AREA_FORMATION


DATABASE_FILE = path.get_full_path(__file__, "data", "profiles.hmm")
RULE_FILE = path.get_full_path(__file__, "cluster_rules", "rules.txt")
SIGNATURE_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")

def _build_ruleset(single_rule: str = "") -> Ruleset:
    categories = {"Synthase-dependent", "Sucrase-dependent", "Monosaccharide-synthesis"}  # contains all categories in the rules that will
                           # be used in the ruleset
    
    signatures = {sig.name: sig for sig in get_signature_profiles(SIGNATURE_FILE)}


    rules = create_rules([RULE_FILE], signature_names=set(signatures),
                         valid_categories=categories)
    if single_rule:
        rules = [rule for rule in rules if rule.name == single_rule]

    return Ruleset(
        rules,
        signatures,
        valid_categories=categories,
        database_file=DATABASE_FILE,
        tool=NAME,
        equivalence_groups=[],
    )


_RULESET = _build_ruleset()


class CustomDetectionResults(DetectionResults):
    """ A container for clusters predicted by rules in this module """
    schema_version = 1

    def __init__(self, record_id: str, rule_results: RuleDetectionResults, restricted_to: str) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.restricted_to = restricted_to

    def to_json(self) -> dict[str, Any]:
        # extend this as necessary, covering the full results so it can be regenerated
        return {
            "record_id": self.record_id,
            "schema_version": self.schema_version,
            "restricted_to": self.restricted_to,
            "rule_results": self.rule_results.to_json(),
        }

    @staticmethod
    def from_json(json: dict[str, Any], record: Record) -> Self:
        # checking the input is valid is a good idea, but is omitted here
        rule_results = RuleDetectionResults.from_json(json["rule_results"], record)
        if rule_results is None:
            raise ValueError("Detection results have changed. No results can be reused")

        return CustomDetectionResults(
            json["record_id"],
            rule_results,
            json["restricted_to"],
        )

    def get_predicted_protoclusters(self) -> list[Protocluster]:
        """ Used by core antiSMASH logic to add protoclusters to the record """
        return self.rule_results.protoclusters


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module, these
        will automatically be included by the usual antiSMASH module handling.
    """
    # starting with the base, supplying the description and prefix for all options
    args = ModuleArgs(
        title="Trivial detection module",  # the name of the argument group in which
                                           # to bundle arguments in the help output
        prefix="det",  # prefix for all commands added below
        enabled_by_default=True,  # whether this module should always run
        basic_help=False,  # whether to show the help for this module in the
                           # basic '--help' or only in the full '--help-showall'
    )
    # add a toggle for this module, specifically to disable it as it is enabled
    # by default above
    args.add_analysis_toggle(
        'disable',   # the commmand line argument itself (the prefix is added automatically)
         dest='disabled',  # the naming of the result in the options object (again prefix is added)
         default=False,
         action='store_true',
         help="Run TFBS finder on all gene clusters."
     )
    # any other options that will take particular values, rather than simple on/off
    args.add_option(
        "single-rule",
        dest="single_rule",
        # for remaining options, see the documentation for the 'argparse' standard library module
        type=str,
        default="",
        help=("Resticts the rules used to the given rule name")
    )

    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return not options.det_disabled


def run_on_record(record: Record, previous_results: Optional[CustomDetectionResults],
                  options: ConfigType) -> CustomDetectionResults:
    """ This is where the analysis itself happens, running over the record and
        and generating results.
    """
    if previous_results:
        return previous_results

    ruleset = _build_ruleset(options.det_single_rule)
    if options.det_single_rule:
        logging.info("detection restricted to: %s", options.det_single_rule)
    results = detect_protoclusters_and_signatures(record, ruleset)
    results.annotate_cds_features()
    return CustomDetectionResults(record.id, results, restricted_to=options.det_single_rule)


def regenerate_previous_results(results: dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[CustomDetectionResults]:
    """ This would normally rebuild any results from a JSON-friendly format, for
        when the '--reuse' option is supplied on the command line.
    """
    # options should be checked here to see if they were changed from the previous results,
    # but this step is omitted in the demo
    return CustomDetectionResults.from_json(results, record)


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures the module data is fully prepared, e.g. HMM profile database is pressed """
    failure_messages = []

    failure_messages.extend(
        ensure_database_pressed(
            DATABASE_FILE,
            return_not_raise=True,  # allows all errors to be described, not just the first
        )
    )
    # any other data used should be prepared here, e.g. sklearn classifiers

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check that all prerequistes are satisfied, e.g. binary dependencies and
        datafiles.
    """
    # for this specific demo module, it will reuse the check from antiSMASH's hmm_detection
    return original_check_prereqs(options)


def check_options(options: ConfigType) -> list[str]:
    """ Check that all options are valid """
    failure_messages = []
    # the one option defined is to restrict the ruleset down to a single rule
    # if that option isn't in the rules, that's an error
    if options.det_single_rule:
        try:
            _RULESET.get_rule_by_name(options.det_single_rule)
        except ValueError:
            failure_messages.append(f"Ruleset '{options.det_single_rule}' does not exist")

    # any other options should also be checked here

    return failure_messages
