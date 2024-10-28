# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This is an arbitrarily tiny analysis module that keeps track of any CDS with a specific string
in the identifier. For more detailed examples, see the template in the antiSMASH
repository's wiki or actual modules within antiSMASH.

See the included customsmash.detection.custom_detection module for more detailed information.

No HTML sections are generated, nothing is added to the record for genbank output.
Only the JSON output will have these results.
"""

from typing import Any, Optional, Self

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

NAME = "custom_analysis"
SHORT_DESCRIPTION = "Keeps track of genes with a specific name fragment"

class NameResults(ModuleResults):
    def __init__(self, record_id: str, target: str, names: str):
        super().__init__(record_id)
        self.target = target
        self.names = names

    def add_to_record(self, record: Record) -> None:
        # normally any results that need to appear in genbank files would
        # be added here, but for now skip it
        return

    def to_json(self) -> dict[str, Any]:
        return {
            "record_id": self.record_id,
            "target": self.target,
            "names": self.names,
        }

    @staticmethod
    def from_json(json: dict[str, Any], record: Record) -> Self:
        # this is a cheat to regenerate the results when the JSON naming
        # matches the variable naming *and* all types are simple python types
        return NameResults(**json)


def get_arguments() -> ModuleArgs:
    args = ModuleArgs("Additional analysis", 'naming', enabled_by_default=True)
    args.add_option("target",
                    dest="target",
                    type=str,
                    default="",
                    help="A string fragment to track within gene names")
    return args


def check_options(_options: ConfigType) -> list[str]:
    # the options are trivial here, so there won't be issues
    return []


def check_prereqs(_options: ConfigType) -> list[str]:
    # no dependencies of any kind
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return not options.minimal  # a pipeline-wide option to disable everything unnecessary


def regenerate_previous_results(previous: dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[NameResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return NameResults.from_json(previous, record)


def run_on_record(record: Record, results: Optional[NameResults], options: ConfigType) -> NameResults:
    """ Run the analysis, unless the previous results apply to the given record """
    if results:
        return results

    if not options.naming_target:
        return NameResults(record.id, "", [])

    found = []
    for cds in record.get_cds_features_within_regions():
        if options.naming_target in cds.get_name():
            found.append(cds)
    return NameResults(record.id, options.target, found)
