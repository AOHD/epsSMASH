# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML generation module using antiSMASH HTML as a base
"""

import argparse
import glob
import logging
import os
import re
import shutil
from typing import Any, Dict, Iterable, List, Optional, Tuple
import warnings
import csv
import importlib.resources
import sass

from antismash.outputs import html
from antismash.outputs.html import copy_template_dir
from antismash.outputs.html.generator import Legend
from antismash.outputs.html import generator

from antismash.common import html_renderer, path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Feature, Record, Region
from antismash.custom_typing import AntismashModule
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.outputs.html.generator import generate_webpage, find_local_antismash_js_path
from antismash.modules import clusterblast

NAME = "html"
SHORT_DESCRIPTION = "HTML output"
print("TEMPLATE_PATH:", generator.TEMPLATE_PATH)
generator.TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "templates")

# Load legends to override the default ones

LEGENDS = [
    Legend(css_class="legend-type-biosynthetic", label="core biosynthetic genes"),
    Legend(css_class="legend-type-biosynthetic-additional", label="additional biosynthetic genes"),
    Legend(css_class="legend-type-transport", label="transport related genes"),
    Legend(css_class="legend-type-regulatory", label="regulatory genes"),
    Legend(css_class="legend-type-degradation", label="degradation genes"),
    Legend(css_class="legend-type-glycosyltransferase", label="glycosyltransferase genes"),
    Legend(css_class="legend-type-modification", label="modification genes"),
    Legend(css_class="legend-type-polymerisation", label="polymerisation genes"),
    Legend(css_class="legend-type-other", label="other genes"),
    Legend(css_class="legend-type-precursor", label="precursor genes")
]

original_convert = html.js.convert_cds_features

# load table with detection profile names matched to function types 
# e.g. a file with "profile name\ttransport" 

def read_tsv_to_dict():
    table = {}
    with importlib.resources.open_text("epssmash.outputs.html", "gene_functions.tsv") as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) == 2:
                key, value = row
                table[key] = value
    return table

def format_dict_as_string(table):
    formatted_string = "#TABLE = {\n"
    for key, value in table.items():
        formatted_string += f'    "{key}": "{value}",\n'
    formatted_string += "}"
    return formatted_string


TABLE = read_tsv_to_dict()

def convert_cds_features(record: Record, features: Iterable[CDSFeature], options: ConfigType,
                         mibig_entries: Dict[str, List[clusterblast.results.MibigEntry]], offset: int = 0,
                         ) -> List[Dict[str, Any]]:
    """ Convert CDSFeatures to JSON """
    original_results = original_convert(record, features, options, mibig_entries) 

    for feature, js in zip(features, original_results):
        
        detection_results = feature.gene_functions.get_by_tool("rule-based-clusters")
        if len(detection_results) == 0:
            continue
        
        else:
            # replace function for the purposes of colours
            # TODO: don't just use the first, figure out how to pick which of many you might want to use
            # this will also override CORE functions
            js["type"] = TABLE.get(detection_results[0].description, "other")
            # the type needs to be present in CSS for both ".svgene-type-<new function type>"
            # and ".legend-type-<function type>"
    
    return original_results


# replace the original with the wrapper
html.js.convert_cds_features = convert_cds_features


def get_arguments() -> ModuleArgs:
    """ Builds the arguments for the HMTL output module """
    # shortcut here to use the antiSMASH HTML arguments, but they could be
    # replaced completely or just added to
    args = html.get_arguments()
    return args


def prepare_data(_logging_only: bool = False) -> List[str]:
    """ Rebuild any dynamically buildable data """
    flavours = ["bacteria"]

    with path.changed_directory(path.get_full_path(__file__, "css")):
        built_files = [os.path.abspath(f"{flavour}.css") for flavour in flavours]

        if path.is_outdated(built_files, glob.glob("*.scss")):
            logging.info("CSS files out of date, rebuilding")

            for flavour in flavours:
                target = f"{flavour}.css"
                source = f"{flavour}.scss"
                assert os.path.exists(source), flavour
                result = sass.compile(filename=source, output_style="compact")
                with open(target, "w", encoding="utf-8") as out:
                    out.write(result)
    return []


def check_prereqs(_options: ConfigType) -> List[str]:
    """ Check prerequisites """
    return prepare_data()


def check_options(_options: ConfigType) -> List[str]:
    """ Check options, but none to check here """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Is the HMTL module enabled (currently always enabled) """
    return options.html_enabled or not options.minimal


def write(records: List[Record], results: List[Dict[str, ModuleResults]],
          options: ConfigType, all_modules: List[AntismashModule]) -> None:
    """ Writes all results to a webpage, where applicable. Writes to options.output_dir

        Arguments:
            records: the list of Records for which results exist
            results: a list of dictionaries containing all module results for records
            options: antismash config object
            all_modules: a list of all modules which might create sections of HTML

        Returns:
            None
    """
    output_dir = options.output_dir

    copy_template_dir(path.get_full_path(__file__, "css"), output_dir, pattern=f"{options.taxon}.css")
    # reuse the antiSMASH default javascript
    # if modifications are required, then provide the javascript files and change the path here
    copy_template_dir(path.get_full_path(html.__file__, "js"), output_dir)
    # if there wasn't an antismash.js in the JS dir, fall back to one in databases
    local_path = os.path.join(output_dir, "js", "antismash.js")
    if not os.path.exists(local_path):
        js_path = find_local_antismash_js_path(options)
        if js_path:
            logging.debug("Results page using antismash.js from local copy: %s", js_path)
            shutil.copy(js_path, os.path.join(output_dir, "js", "antismash.js"))
    # and if it's still not there, that's fine, it'll use a web-accessible URL
    if not os.path.exists(local_path):
        logging.debug("Results page using antismash.js from remote host")

    # copy non-antismash specific images from antismash proper
    copy_template_dir(path.get_full_path(html.__file__, "images"), output_dir)
    # and then all the replacements and/or additions
    copy_template_dir(path.get_full_path(__file__, "images"), output_dir, keep_existing_content=True)

    with open(os.path.join(options.output_dir, "index.html"), "w", encoding="utf-8") as result_file:
        content = generate_webpage(records, results, options, all_modules, legends=LEGENDS)
        # strip all leading whitespace and blank lines, as they're meaningless to HTML
        content = re.sub("^( *|$)", "", content, flags=re.M)
        result_file.write(content)
