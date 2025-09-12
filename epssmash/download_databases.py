# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to download epsSMASH data files not shipped with the code."""

import argparse
import gzip
import hashlib
import lzma
import os
import pathlib
import sys
import tarfile
from typing import Any, Type
from urllib import error as urlerror
from urllib import request

import antismash
from antismash.common.hmmer import ensure_database_pressed
from antismash.common.html_renderer import (
    get_antismash_js_version,
    get_antismash_js_url,
)
from antismash.common import json

PFAM_LATEST_VERSION = "35.0"
PFAM_LATEST_URL = f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{PFAM_LATEST_VERSION}/Pfam-A.hmm.gz"
PFAM_LATEST_ARCHIVE_CHECKSUM = "48ec2d1123c84046b00279eae1fb3d5be1b578e6221453f329d16954c89d0d35"
PFAM_LATEST_CHECKSUM = "8d3e2ffa785f91ee0e24a3994d2dcfff6f382e3cf663784a47688e7d95297fee"

CLUSTERBLAST_URL = "https://dl.secondarymetabolites.org/releases/epssmash/clusterblast/clusterblast_1.0.tar.xz"
CLUSTERBLAST_ARCHIVE_CHECKSUM = "fc5c8f13292c545fca50a0dcd0b806f83d2af88030b8e487991e6a121892cd51"
CLUSTERBLAST_FASTA_CHECKSUM = "b68d3c47cf191e282a16fb6e386442fe"

LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))

CHUNK = 128 * 1024


class DownloadError(RuntimeError):
    """Exception to throw when downloads fail."""

    pass  # pylint: disable=unnecessary-pass


def get_remote_filesize(url: str) -> int:
    """Get the file size of the remote file."""
    try:
        with request.urlopen(request.Request(url, method="HEAD")) as usock:
            dbfilesize = usock.info().get("Content-Length", "0")
    except urlerror.URLError:
        dbfilesize = "0"

    dbfilesize = int(dbfilesize)  # db file size in bytes
    return dbfilesize


def get_free_space(folder: str) -> int:
    """Return folder/drive free space (in bytes)."""
    return os.statvfs(folder).f_bfree * os.statvfs(folder).f_frsize


def check_diskspace(file_url: str) -> None:
    """Check if sufficient disk space is available."""
    dbfilesize = get_remote_filesize(file_url)
    free_space = get_free_space(".")
    if free_space < dbfilesize:
        raise DownloadError(
            "ERROR: Insufficient disk space available (required: {dbfilesize}, free: {free_space})."
        )


def download_file(url: str, filename: str) -> str:
    """Download a file."""
    try:
        req = request.urlopen(url)  # pylint: disable=consider-using-with
    except urlerror.URLError:
        raise DownloadError("ERROR: File not found on server.\nPlease check your internet connection.")

    # use 1 because we want to divide by the expected size, can't use 0
    expected_size = int(req.info().get("Content-Length", "1"))

    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    overall = 0
    with open(filename, "wb") as handle:
        while True:
            try:
                chunk = req.read(CHUNK)
                if not chunk:
                    print("")
                    break
                overall += len(chunk)
                print(
                    f"\rDownloading {basename}: {(overall / expected_size) * 100:5.2f}% downloaded.",
                    end="",
                )
                handle.write(chunk)
            except IOError:
                raise DownloadError("ERROR: Download interrupted.")
    return filename


def checksum(filename: str, chunksize: int = 2 ** 20) -> str:
    """Get the SHA256 checksum of a file."""
    sha = hashlib.sha256()
    with open(filename, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunksize), b""):
            sha.update(chunk)

    return sha.hexdigest()


def unzip_file(filename: str, decompressor: Any, error_type: Type[Exception]) -> str:
    """Decompress a compressed file."""
    newfilename, _ = os.path.splitext(filename)
    basename = os.path.basename(filename)
    try:
        zipfile = decompressor.open(filename, "rb")
        chunksize = 128 * 1024
        with open(newfilename, "wb") as handle:
            while True:
                try:
                    chunk = zipfile.read(chunksize)
                    if not chunk:
                        break
                    handle.write(chunk)
                except IOError:
                    raise DownloadError("ERROR: Unzipping interrupted.")
    except error_type:
        raise RuntimeError(
            f"Error extracting {basename}. Please try to extract it manually."
        )
    print(f"Extraction of {basename} finished successfully.")
    return newfilename

def untar_file(filename: str) -> None:
    """Extract a TAR/GZ file."""
    basename = filename.rpartition(os.sep)[2]
    try:
        # Remove any version pattern like _X.X or _X.X.X from the directory name
        import re
        base_name = os.path.basename(filename)
        # Remove file extensions
        dir_name = re.sub(r'\.(tar\.xz|tar\.gz|tar)$', '', base_name)
        # Remove version patterns like _1.0, _2.1.3, etc.
        dir_name = re.sub(r'_\d+(\.\d+)*$', '', dir_name)
        extract_dir = os.path.join(os.path.dirname(filename), dir_name)
        
        with tarfile.open(filename) as tar:
            tar.extractall(path=extract_dir)
    except tarfile.ReadError:
        print(
            f"ERROR: Error extracting {basename}. Please try to extract it manually."
        )
        return
    print(f"Extraction of {basename} finished successfully.")

def delete_file(filename: str) -> None:
    """Delete a file."""
    try:
        os.remove(filename)
    except OSError:
        pass


def present_and_checksum_matches(filename: str, sha256sum: str) -> bool:
    """Check if a file is present and the checksum matches."""
    if os.path.exists(filename):
        print(f"Creating checksum of {os.path.basename(filename)}")
        csum = checksum(filename)
        if csum == sha256sum:
            return True
    return False


def download_if_not_present(url: str, filename: str, sha256sum: str) -> None:
    """Download a file if it's not present or checksum doesn't match."""
    # If we are missing the archive file, go and download
    if not present_and_checksum_matches(filename, sha256sum):
        download_file(url, filename)

    print(f"Creating checksum of {os.path.basename(filename)}")
    csum = checksum(filename)
    if csum != sha256sum:
        raise DownloadError(
            f"Error downloading {filename}, sha256sum mismatch. Expected {sha256sum}, got {csum}."
        )


def download_antismash_js(db_dir: str) -> None:
    """ Downloads the latest relevant version of the antiSMASH javascript """
    version = get_antismash_js_version()
    url = get_antismash_js_url()
    download_file(url, os.path.join(db_dir, "as-js", version, "antismash.js"))


def download_pfam(db_dir: str, url: str, version: str, archive_checksum: str, db_checksum: str) -> None:
    """Download and compile the PFAM database."""
    archive_filename = os.path.join(db_dir, "pfam", version, "Pfam-A.hmm.gz")
    db_filename = os.path.splitext(archive_filename)[0]

    if present_and_checksum_matches(db_filename, db_checksum):
        print("PFAM file present and ok for version", version)
        return

    print("Downloading PFAM version", version)
    check_diskspace(url)
    download_if_not_present(url, archive_filename, archive_checksum)
    filename = unzip_file(archive_filename, gzip, gzip.zlib.error)  # type: ignore
    ensure_database_pressed(filename)
    delete_file(filename + ".gz")


def download_clusterblast(db_dir: str) -> None:
    """Download the clusterblast database."""
    archive_filename = os.path.join(db_dir, CLUSTERBLAST_URL.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "clusterblast", "proteins.fasta")

    if present_and_checksum_matches(fasta_filename, CLUSTERBLAST_FASTA_CHECKSUM):
        print("ClusterBlast fasta file present and checked")
        return

    print("Downloading ClusterBlast database.")
    check_diskspace(CLUSTERBLAST_URL)
    download_if_not_present(CLUSTERBLAST_URL, archive_filename, CLUSTERBLAST_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download(args: argparse.Namespace) -> bool:
    """Download all the large external databases needed, along with non-python
       components.
    """
    download_antismash_js(args.database_dir)

    # grab the latest pfam
    download_pfam(
        args.database_dir,
        PFAM_LATEST_URL,
        PFAM_LATEST_VERSION,
        PFAM_LATEST_ARCHIVE_CHECKSUM,
        PFAM_LATEST_CHECKSUM,
    )

    download_clusterblast(args.database_dir)
    
    return True


def _main() -> None:
    """ Downloads, decompresses, and compiles large databases. Also ensures
        antiSMASH's module data is prepared.
    """
    # Small dance to grab the antiSMASH config for the database dir.
    # All the modules are required to parse the config file,
    # and any executable paths defined should be kept.
    all_modules = antismash.get_detection_modules() + antismash.get_analysis_modules()
    config = antismash.config.build_config(args=[], parser=None, isolated=False, modules=all_modules)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database-dir",
        default=config.database_dir,
        metavar="DIR",
        help="Base directory for the antiSMASH databases (default: %(default)s).",
    )

    args = parser.parse_args()
    antismash.config.update_config({"database_dir": args.database_dir})
    if not download(args):
        print("Errors occurred while downloading, aborted.")
        sys.exit(1)
    try:
        print("Pre-building all databases...")
        antismash.main.prepare_module_data()
        print("done.")
    except Exception as err:  # pylint: disable=broad-except
        print("Error encountered while preparing module data:", str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _main()