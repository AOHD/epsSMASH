# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""

"""

import sys

import antismash
from antismash.__main__ import main as antismash_main
# import some modules from antiSMASH that won't be overridden
# in particula some generic module
from antismash.detection import cluster_hmmer, genefunctions
from antismash.modules import pfam2go

# import those modules that are defined only within this custom pipeline,
# the structure and naming here are just arbitrary examples and can be changed
# to suit
from customsmash.detection import custom_detection
from customsmash.modules import custom_analysis
from customsmash.outputs import html

# replace the normal antismash modules with any combination of antiSMASH modules
# and/or custom modules
antismash.main.replace_analysis_modules([custom_analysis, pfam2go])
antismash.main.replace_detection_modules([cluster_hmmer, custom_detection, genefunctions])
antismash.main.replace_html_module(html)

# override search path for any user config file
# a good naming convention is <name><major version> to avoid conflicts between version options
antismash.config.set_user_config_file("~/.customsmash1.cfg")  # if it doesn't exist, it's ignored
# override default config file path
antismash.config.set_alternate_defaults_file(antismash.common.path.get_full_path(__file__, "config", "default.cfg"))

__version__ = "1.0"


def main(args: list[str]) -> int:
    # Any custom processing before or after running the pipeline can be added here

    # e.g. some extra console printing
    print("CustomSMASH: a pipeline that finds genes with identifiers that begin with 'A'")
    # The simplest version of running is to use antiSMASH's own main which will handle
    # all argument parsing. Alternatively, use antismash.main.run_antismash()
    # with the sequence file and the config object built via
    # antismash.config.build_config(). In either case, any command line arguments
    # defined by included modules will be handled automatically.

    # using main won't give much information about what went wrong, if it went wrong
    # overriding the name and version here will propagate to all those places where the name "antiSMASH"
    # would be used
    result_code = antismash_main(args, branding_override="CustomSMASH", version_override=__version__)

    # but any following processing can still take place, e.g.
    print("CustomSMASH has finished")

    return result_code


if __name__ == "__main__":
    sys.exit(main(sys.argv))
