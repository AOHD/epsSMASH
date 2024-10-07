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
import os
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
from antismash.common.hmm_rule_parser.structures import Multipliers
from antismash.config.args import ModuleArgs, SplitCommaAction

NAME = "customsmash_detection"
SHORT_DESCRIPTION = "some kind of protocluster detection"
# the detection stage defines when the module is run in the detection process
DETECTION_STAGE = DetectionStage.AREA_FORMATION


HMM_FILE = path.get_full_path(__file__, "data", "bgc_seeds.hmm")
#RULE_FILE = path.get_full_path(__file__, "cluster_rules", "rules.txt")
SIGNATURE_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")

_STRICTNESS_LEVELS = ["strict", "relaxed", "loose"]


_RULESETS: dict[tuple[str, tuple[str, ...], tuple[str, ...], Multipliers], Ruleset] = {}

def _get_rule_files_for_strictness(strictness: str) -> list[str]:
    """ Returns a list of appropriate rule files for the given strictness level """
    assert strictness in _STRICTNESS_LEVELS, strictness
    files = []
    for level in _STRICTNESS_LEVELS[:_STRICTNESS_LEVELS.index(strictness) + 1]:
        files.append(path.get_full_path(__file__, "cluster_rules", f"{level}.txt"))
    return files


def _build_ruleset(options: ConfigType) -> Ruleset:
    

    strictness = options.hmmdetection_strictness
    name_subset = set(options.hmmdetection_limit_to_rules)
    category_subset = set(options.hmmdetection_limit_to_categories)
    
    # the cache key needs to be immutable
    key = (strictness, tuple(name_subset), tuple(category_subset))
    
    categories = {"Synthase-dependent", "Sucrase-dependent", "Monosaccharide-synthesis", "Wzy-dependent"}  # contains all categories in the rules that will
                           # be used in the ruleset
    
    signatures = {sig.name: sig for sig in get_signature_profiles(SIGNATURE_FILE)}

     # return any existing ruleset
    ruleset = _RULESETS.get(key)
    if ruleset:
        return ruleset

    # otherwise make a default ruleset for the strictness
    ruleset = Ruleset.from_files(
        signature_file = SIGNATURE_FILE, 
        seeds = HMM_FILE,
        rule_files = _get_rule_files_for_strictness(strictness),
        categories = categories,
        filter_file = os.devnull,
        tool = "rule-based-clusters")

    # limit the rules used, if relevant
    rules: Iterable[rule_parser.DetectionRule] = ruleset.rules
    if name_subset:
        rules = filter(lambda rule: rule.name in name_subset, rules)
    if category_subset:
        rules = filter(lambda rule: rule.category in category_subset, rules)

    ruleset = ruleset.copy_with_replacements(rules=list(rules))

    # update the cache
    _RULESETS[key] = ruleset

    return ruleset


class CustomDetectionResults(DetectionResults):
    """ A container for clusters predicted by rules in this module """
    schema_version = 1

    def __init__(self, record_id: str, rule_results: RuleDetectionResults, restricted_to: list[str],
                 strictness: str) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.restricted_to = restricted_to
        if strictness not in _STRICTNESS_LEVELS:
            raise ValueError(f"unknown strictness level: {strictness}")
        self.strictness = strictness

    def to_json(self) -> dict[str, Any]:
        # extend this as necessary, covering the full results so it can be regenerated
        return {
            "record_id": self.record_id,
            "schema_version": self.schema_version,
            "restricted_to": self.restricted_to,
            "rule_results": self.rule_results.to_json(),
            "strictness": self.strictness
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
            json["strictness"]
        )

    def get_predicted_protoclusters(self) -> list[Protocluster]:
        """ Used by core antiSMASH logic to add protoclusters to the record """
        return self.rule_results.protoclusters


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs('HMM detection options', 'hmmdetection')
    # add a toggle for this module, specifically to disable it as it is enabled
    # by default above
    args.add_analysis_toggle(
        'disable',   # the commmand line argument itself (the prefix is added automatically)
         dest='disabled',  # the naming of the result in the options object (again prefix is added)
         default=False,
         action='store_true',
         help="Run cluster detection."
     )
    args.add_option('strictness',
                    dest='strictness',
                    type=str,
                    choices=["strict", "relaxed", "loose"],
                    default="loose",
                    help=("Defines which level of strictness to use for "
                          "HMM-based cluster detection, (default: %(default)s)."))
    args.add_option("limit-to-rule-names",
                    dest="limit_to_rules",
                    metavar="RULE1[,RULE2,...]",
                    action=SplitCommaAction,
                    default=[],
                    help="Restrict detection to the named rules (default: no limits).")
    args.add_option("limit-to-rule-categories",
                    dest="limit_to_categories",
                    metavar="CATEGORY1[,CATEGORY2,...]",
                    action=SplitCommaAction,
                    default=[],
                    help="Restrict detection to the given rules (default: no limits).")

    return args

def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return not options.hmmdetection_disabled


def run_on_record(record: Record, previous_results: Optional[CustomDetectionResults],
                  options: ConfigType) -> CustomDetectionResults:
    """ This is where the analysis itself happens, running over the record and
        and generating results.
    """
    if previous_results:
        return previous_results
    
    strictness = options.hmmdetection_strictness
    logging.info("HMM detection using strictness: %s", strictness)


    ruleset = _build_ruleset(options)
    if options.hmmdetection_limit_to_rules:
        logging.info("detection restricted to: %s", options.hmmdetection_limit_to_rules)
    
    if options.hmmdetection_strictness:
        logging.info("detection strictness: %s", options.hmmdetection_strictness)
    
    results = detect_protoclusters_and_signatures(record, ruleset)
    results.annotate_cds_features()
    return CustomDetectionResults(record.id, results, restricted_to=options.hmmdetection_limit_to_rules, strictness=strictness)


def regenerate_previous_results(results: dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[CustomDetectionResults]:
    """ This would normally rebuild any results from a JSON-friendly format, for
        when the '--reuse' option is supplied on the command line.
    """
    # options should be checked here to see if they were changed from the previous results,
    # but this step is omitted in the demo
    return CustomDetectionResults.from_json(results, record)


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages = []

    # Check that hmmdetails.txt is readable and well-formatted
    try:
        profiles = get_signature_profiles(SIGNATURE_FILE)
    except ValueError as err:
        if not logging_only:
            raise
        return [str(err)]

    # the path to the markov model
    seeds_hmm = path.get_full_path(__file__, 'data', 'bgc_seeds.hmm')
    hmm_files = [os.path.join("data", "individual_hmms", sig.hmm_file) for sig in profiles]
    # include the listing, since tools like wget will keep modified timestamps on the HMMs
    description_file = path.get_full_path(__file__, 'data', 'hmmdetails.txt')
    outdated = False
    if not path.locate_file(seeds_hmm):
        logging.debug("%s: %s doesn't exist, regenerating", NAME, seeds_hmm)
        outdated = True
    else:
        seeds_timestamp = os.path.getmtime(seeds_hmm)
        for component in hmm_files + [description_file]:
            if os.path.getmtime(component) > seeds_timestamp:
                logging.debug("%s out of date, regenerating", seeds_hmm)
                outdated = True
                break

    # regenerate if missing or out of date
    if outdated:
        # try to generate file from all specified profiles in hmmdetails
        try:
            with open(seeds_hmm, "w", encoding="utf-8") as all_hmms_handle:
                for hmm_file in hmm_files:
                    with open(path.get_full_path(__file__, hmm_file), "r", encoding="utf-8") as handle:
                        all_hmms_handle.write(handle.read())
        except OSError:
            if not logging_only:
                raise
            failure_messages.append(f"Failed to generate file {seeds_hmm!r}")

    # if regeneration failed, don't try to run hmmpress
    if failure_messages:
        return failure_messages

    failure_messages.extend(ensure_database_pressed(seeds_hmm, return_not_raise=logging_only))

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check that all prerequistes are satisfied, e.g. binary dependencies and
        datafiles.
    """
    # for this specific demo module, it will reuse the check from antiSMASH's hmm_detection
    return prepare_data() + original_check_prereqs(options)


def check_options(options: ConfigType) -> list[str]:
    """ Check that all options are valid """
    failure_messages = []
    # the one option defined is to restrict the ruleset down to a single rule
    # if that option isn't in the rules, that's an error
    if options.hmmdetection_limit_to_rules:
        try:
            _RULESET.get_rule_by_name(options.hmmdetection_limit_to_rules)
        except ValueError:
            failure_messages.append(f"Ruleset '{options.hmmdetection_limit_to_rules}' does not exist")
    if options.hmmdetection_strictness not in _STRICTNESS_LEVELS:
        issues.append(f"Unknown strictness level: {options.strictness}")

    # any other options should also be checked here

    return failure_messages
