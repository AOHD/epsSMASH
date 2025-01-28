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

import sass

from antismash.outputs import html
from antismash.outputs.html import copy_template_dir
from antismash.outputs.html.generator import Legend

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

#TABLE = {
#    "eg_acetan": "degradation",
#    "gumE_acetan": "polymerisation",
#    "gumF_acetan": "modification",
#    "gumK_acetan": "glycosyltransferase",
#    "aceR_acetan": "glycosyltransferase",
#    "aceQ_acetan": "glycosyltransferase",
#    "aceP_acetan": "glycosyltransferase",
#    "aceM_acetan": "precursor",
#    "aceF_acetan": "precursor",
#    "aceA_acetan": "glycosyltransferase",
#    "aceB_acetan": "glycosyltransferase",
#    "aceD_acetan": "transport",
#    "aceC_acetan": "glycosyltransferase",
#    "aceE_acetan": "transport",
#    "aceG_acetan": "transport",
#    "aceH_acetan": "transport",
#    "aceI_acetan": "modification",
#    "epsA_lactobacillus": "regulatory",
#    "epsB_lactobacillus": "transport",
#    "epsC_lactobacillus": "transport",
#    "epsD_lactobacillus": "regulatory",
#    "epsE_lactobacillus": "glycosyltransferase",
#    "1026_lactobacillus": "glycosyltransferase",
#    "1027_lactobacillus": "glycosyltransferase",
#    "1028_lactobacillus": "glycosyltransferase",
#    "1029_lactobacillus": "glycosyltransferase",
#    "1030_lactobacillus": "glycosyltransferase",
#    "1031_lactobacillus": "polymerisation",
#    "1032_lactobacillus": "modification",
#    "1033_lactobacillus": "glycosyltransferase",
#    "1034_lactobacillus": "transport",
#    "epsA_pnag": "regulatory",
#    "epsB_pnag": "regulatory",
#    "epsC_pnag": "glycosyltransferase",
#    "epsD_pnag": "glycosyltransferase",
#    "epsE_pnag": "glycosyltransferase",
#    "epsF_pnag": "glycosyltransferase",
#    "epsG_pnag": "biosynthetic-additional",
#    "epsH_pnag": "glycosyltransferase",
#    "epsI_pnag": "modification",
#    "epsJ_pnag": "glycosyltransferase",
#    "epsK_pnag": "transport",
#    "epsL_pnag": "precursor",
#    "epsM_pnag": "precursor",
#    "epsN_pnag": "precursor",
#    "epsO_pnag": "glycosyltransferase",
#    "pslA_psl": "glycosyltransferase",
#    "pslB_psl": "precursor",
#    "pslC_psl": "glycosyltransferase",
#    "pslD_psl": "transport",
#    "pslE_psl": "transport",
#    "pslF_psl": "glycosyltransferase",
#    "pslG_psl": "degradation",
#    "pslH_psl": "glycosyltransferase",
#    "pslI_psl": "glycosyltransferase",
#    "pslJ_psl": "polymerisation",
#    "pslK_psl": "transport",
#    "pslL_psl": "modification",
#    "exsA_succinoglycan": "transport",
#    "exsB_succinoglycan": "regulatory",
#    "exsC_succinoglycan": "biosynthetic-additional",
#    "exsD_succinoglycan": "biosynthetic-additional",
#    "exsE_succinoglycan": "transport",
#    "exsF_succinoglycan": "regulatory",
#    "exsG_succinoglycan": "regulatory",
#    "exsH_succinoglycan": "degradation",
#    "exsI_succinoglycan": "regulatory",
#    "exoA_succinoglycan": "glycosyltransferase",
#    "exoB_succinoglycan": "precursor",
#    "exoF_succinoglycan": "precursor",
#    "exoH_succinoglycan": "modification",
#    "exoI_succinoglycan": "degradation",
#    "exoK_succinoglycan": "degradation",
#    "exoL_succinoglycan": "glycosyltransferase",
#    "exoM_succinoglycan": "glycosyltransferase",
#    "exoN_succinoglycan": "precursor",
#    "exoO_succinoglycan": "glycosyltransferase",
#    "exoP_succinoglycan": "transport",
#    "exoQ_succinoglycan": "polymerisation",
#    "exoT_succinoglycan": "transport",
#    "exoU_succinoglycan": "glycosyltransferase",
#    "exoV_succinoglycan": "modification",
#    "exoW_succinoglycan": "glycosyltransferase",
#    "exoX_succinoglycan": "regulatory",
#    "exoY_succinoglycan": "precursor",
#    "exoZ_succinoglycan": "modification",
#    "wgeH_galactoglucan": "biosynthetic-additional",
#    "wgeG_galactoglucan": "glycosyltransferase",
#    "wgeF_galactoglucan": "biosynthetic-additional",
#    "wgeE_galactoglucan": "biosynthetic-additional",
#    "wgeD_galactoglucan": "glycosyltransferase",
#    "wgeC_galactoglucan": "modification",
#    "wgeB_galactoglucan": "glycosyltransferase",
#    "wgeA_galactoglucan": "modification",
#    "wgdB_galactoglucan": "polymerisation",
#    "wgdA_galactoglucan": "transport",
#    "wggR_galactoglucan": "regulatory",
#    "wgcA_galactoglucan": "glycosyltransferase",
#    "wgaA_galactoglucan": "biosynthetic-additional",
#    "wgaB_galactoglucan": "glycosyltransferase",
#    "wgaC_galactoglucan": "glycosyltransferase",
#    "wgaD_galactoglucan": "biosynthetic-additional",
#    "wgaE_galactoglucan": "biosynthetic-additional",
#    "wgaF_galactoglucan": "biosynthetic-additional",
#    "wgaG_galactoglucan": "precursor",
#    "wgaH_galactoglucan": "precursor",
#    "wgaI_galactoglucan": "precursor",
#    "wgaJ_galactoglucan": "precursor",
#    "gumB_xanthan": "transport",
#    "gumC_xanthan": "transport",
#    "gumD_xanthan": "glycosyltransferase",
#    "gumE_xanthan": "polymerisation",
#    "gumF_xanthan": "modification",
#    "gumG_xanthan": "modification",
#    "gumH_xanthan": "glycosyltransferase",
#    "gumI_xanthan": "glycosyltransferase",
#    "gumJ_xanthan": "transport",
#    "gumK_xanthan": "glycosyltransferase",
#    "gumL_xanthan": "modification",
#    "gumM_xanthan": "glycosyltransferase",
#    "wza_colanic": "transport",
#    "wzb_colanic": "transport",
#    "wzc_colanic": "polymerisation",
#    "wcaA_colanic": "glycosyltransferase",
#    "wcaB_colanic": "modification",
#    "wcaC_colanic": "glycosyltransferase",
#    "wcaD_colanic": "polymerisation",
#    "wcaE_colanic": "glycosyltransferase",
#    "wcaF_colanic": "modification",
#    "gmd_colanic": "precursor",
#    "wcaG_colanic": "precursor",
#    "wcaH_colanic": "biosynthetic-additional",
#    "wcaI_colanic": "glycosyltransferase",
#    "manC_colanic": "precursor",
#    "manB_colanic": "precursor",
#    "wcaJ_colanic": "glycosyltransferase",
#    "wzx_colanic": "transport",
#    "wcaK_colanic": "biosynthetic-additional",
#    "wcaL_colanic": "glycosyltransferase",
#    "amsG_amylovoran": "glycosyltransferase",
#    "amsH_amylovoran": "transport",
#    "amsI_amylovoran": "regulatory",
#    "amsA_amylovoran": "transport",
#    "amsB_amylovoran": "glycosyltransferase",
#    "amsC_amylovoran": "polymerisation",
#    "amsD_amylovoran": "glycosyltransferase",
#    "amsE_amylovoran": "glycosyltransferase",
#    "amsF_amylovoran": "degradation",
#    "amsJ_amylovoran": "modification",
#    "amsK_amylovoran": "glycosyltransferase",
#    "amsL_amylovoran": "transport",
#    "galF_amylovoran": "precursor",
#    "galE_amylovoran": "precursor",
#    "spnS_sphingans": "transport",
#    "spnG_sphingans": "polymerisation",
#    "spnR_sphingans": "transport",
#    "spnQ_sphingans": "glycosyltransferase",
#    "spnI_sphingans": "biosynthetic-additional",
#    "spnK_sphingans": "glycosyltransferase",
#    "spnL_sphingans": "glycosyltransferase",
#    "spnJ_sphingans": "biosynthetic-additional",
#    "spnF_sphingans": "biosynthetic-additional",
#    "spnD_sphingans": "transport",
#    "spnC_sphingans": "transport",
#    "spnE_sphingans": "transport",
#    "spnM_sphingans": "transport",
#    "spnN_sphingans": "transport",
#    "atrD_sphingans": "biosynthetic-additional",
#    "atrB_sphingans": "biosynthetic-additional",
#    "spnB_sphingans": "glycosyltransferase",
#    "rmlA_sphingans": "precursor",
#    "rmlC_sphingans": "precursor",
#    "rmlB_sphingans": "precursor",
#    "rmlD_sphingans": "precursor",
#    "urf314_sphingans": "modification",
#    "urf31_sphingans": "modification",
#    "urf34_sphingans": "modification",
#    "wzc_emulsan": "transport",
#    "wzb_emulsan": "regulatory",
#    "wza_emulsan": "transport",
#    "weeA_emulsan": "precursor",
#    "weeB_emulsan": "precursor",
#    "weeC_emulsan": "modification",
#    "wzx_emulsan": "transport",
#    "wzy_emulsan": "polymerisation",
#    "weeD_emulsan": "glycosyltransferase",
#    "weeE_emulsan": "biosynthetic-additional",
#    "weeF_emulsan": "biosynthetic-additional",
#    "weeG_emulsan": "glycosyltransferase",
#    "weeH_emulsan": "glycosyltransferase",
#    "weeI_emulsan": "modification",
#    "weeJ_emulsan": "precursor",
#    "weeK_emulsan": "precursor",
#    "galU_emulsan": "biosynthetic-additional",
#    "ugd_emulsan": "biosynthetic-additional",
#    "pgi_emulsan": "biosynthetic-additional",
#    "galE_emulsan": "precursor",
#    "asnB_zooglan": "regulatory",
#    "zooTRP_zooglan": "biosynthetic-additional",
#    "wzxC_zooglan": "transport",
#    "zooGT8_zooglan": "glycosyltransferase",
#    "zooGT7_zooglan": "glycosyltransferase",
#    "zooGH_zooglan": "degradation",
#    "zooM1_zooglan": "modification",
#    "zooGT6_zooglan": "glycosyltransferase",
#    "zooGT5_zooglan": "glycosyltransferase",
#    "epsH_zooglan": "polymerisation",
#    "wzy_zooglan": "polymerisation",
#    "capK_zooglan": "biosynthetic-additional",
#    "zooGT4_zooglan": "glycosyltransferase",
#    "zooGT3_zooglan": "glycosyltransferase",
#    "asnH_zooglan": "regulatory",
#    "zooGT2_zooglan": "glycosyltransferase",
#    "zooGT1_zooglan": "glycosyltransferase",
#    "zooM2_zooglan": "modification",
#    "zooM3_zooglan": "modification",
#    "zooSA_zooglan": "transport",
#    "zooP_zooglan": "transport",
#    "wzc_zooglan": "transport",
#    "etk_zooglan": "modification",
#    "wza_zooglan": "transport",
#    "wzi_zooglan": "modification",
#    "lolD_zooglan": "transport",
#    "algD": "precursor",
#    "alg8": "polymerisation",
#    "alg44": "regulatory",
#    "algK": "transport",
#    "algE": "transport",
#    "algG": "modification",
#    "algX": "modification",
#    "algL": "degradation",
#    "algI": "modification",
#    "algJ": "modification",
#    "algF": "modification",
#    "algA": "precursor",
#    "Glyco_hydro_68": "polymerisation",
#    "GH70": "polymerisation",
#    "ccsA": "polymerisation",
#    "ccsB": "polymerisation",
#    "ccsZ": "degradation",
#    "ccsI": "modification",
#    "ccsH": "modification",
#    "BcsA": "polymerisation",
#    "BcsB": "polymerisation",
#    "BcsC": "transport",
#    "BcsD": "modification",
#    "BcsE": "regulatory",
#    "BcsF": "regulatory",
#    "BcsG": "modification",
#    "BcsQ": "regulatory",
#    "BcsR": "regulatory",
#    "BcsZ": "degradation",
#    "BcsH": "regulatory",
#    "bglX": "degradation",
#    "BcsO": "biosynthetic-additional",
#    "BcsP": "biosynthetic-additional",
#    "BcsK": "transport",
#    "BcsS": "biosynthetic-additional",
#    "algF_cel": "modification",
#    "algI_cel": "modification",
#    "algJ_cel": "modification",
#    "algX_cel": "modification",
#    "crdA": "transport",
#    "crdS": "polymerisation",
#    "crdC": "transport",
#    "XcsA": "polymerisation",
#    "XcsB": "transport",
#    "XcsC": "degradation",
#    "hasA": "polymerisation",
#    "hasB": "precursor",
#    "hasC": "precursor",
#    "pelD_pos": "regulatory",
#    "pelE_pos": "transport",
#    "pelA_pos": "degradation",
#    "pelF_pos": "polymerisation",
#    "pelG_pos": "modification",
#    "pelA": "degradation",
#    "pelB": "transport",
#    "pelC": "transport",
#    "pelD": "regulatory",
#    "pelE": "transport",
#    "pelF": "polymerisation",
#    "pelG": "modification",
#    "icaA": "polymerisation",
#    "icaB": "modification",
#    "icaC": "polymerisation",
#    "icaD": "regulatory",
#    "pgaA": "transport",
#    "pgaB": "modification",
#    "pgaC": "polymerisation",
#    "pgaD": "regulatory"
#}

def read_tsv_to_dict(tsv_file_path):
    table = {}
    with open(tsv_file_path, mode='r') as file:
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


tsv_file_path = os.path.join(os.path.dirname(__file__), "gene_functions.tsv")
TABLE = read_tsv_to_dict(tsv_file_path)

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
