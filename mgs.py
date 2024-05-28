import json
import urllib.request
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from enum import Enum
from typing import NewType, Optional
from collections import defaultdict
import os
import csv

from pydantic import BaseModel

from pathogen_properties import TaxID
from tree import Tree

BIOPROJECTS_DIR = "bioprojects"


BioProject = NewType("BioProject", str)
Sample = NewType("Sample", str)


target_bioprojects = {
    "crits_christoph": [BioProject("CC-PRJNA661613")],
    "rothman": [BioProject("Rothman-PRJNA729801")],
    "spurbeck": [BioProject("Spurbeck-PRJNA924011")],
    "brinch": [BioProject("Brinch-PRJEB13832"), BioProject("Brinch-PRJEB34633")],
}


class Enrichment(Enum):
    VIRAL = "viral"
    PANEL = "panel"


class SampleAttributes(BaseModel):
    country: str
    state: Optional[str] = None
    county: Optional[str] = None
    location: Optional[str] = None
    fine_location: Optional[str] = None
    # Fixme: Not all the dates are real dates
    date: date | str
    reads: Optional[int] = None
    enrichment: Optional[Enrichment] = None
    method: Optional[str] = None


def european_to_iso(date):
    dd,mm,yyyy = date.split("/")
    return "%s-%s-%s"%(yyyy,mm,dd)

def parse_metadata(record, paper):
    if paper == "rothman":
        sample,library,date,location,enrichment,sample_alias,dataset,bioproject = record
        wtp = sample_alias.split("_")[0] 
        if wtp == "JW":
            # Rothman confirmed over email that JW = JWPCP.
            wtp = "JWPCP"


        return sample, SampleAttributes(
            country = "United States",
            date=date,
            state="California",
            location="Los Angeles",
            county={
                # Hyperion
                "HTP": "Los Angeles County",
                # San Jose Creek
                "SJ": "Los Angeles County",
                # Joint Water Pollution Control Plant
                "JWPCP": "Los Angeles County",
                # Orange County
                "OC": "Orange County",
                # Point Loma
                "PL": "San Diego County",
                # South Bay
                "SB": "San Diego County",
                # North City
                "NC": "San Diego County",
                # Escondido Hale Avenue Resource Recovery Facility
                "ESC": "San Diego County",
            }[wtp],
            fine_location=wtp,
            enrichment="panel" if enrichment == "1" else "viral",
        )
    elif paper == "crits_christoph":
        library,sample,location,date,method,enrichment,sample_alias,dataset,bioproject = record
        return sample, SampleAttributes(
            date=european_to_iso(date),
            country="United States",
            state="California",
            location="San Francisco",
            county={
                "Berkeley": "Alameda County",
                "Marin": "Marin County",
                "Oakland": "Alameda County",
                "SF": "San Francisco County",
            }[location],
            fine_location=location,
            method=method,
            enrichment="panel" if enrichment == "enriched" else "viral",
        )
    elif paper == "spurbeck":
        library,sample,group,date,instrument_model,sample_alias,bioproject,dataset = record
        return sample, SampleAttributes(
            date=european_to_iso(date),
            country="United States",
            state="Ohio",
            location="Ohio",
            # https://github.com/naobservatory/mgs-pipeline/issues/9
            county={
                "A": "Summit County",
                "B": "Trumbull County",
                "C": "Lucas County",
                "D": "Lawrence County",
                "E": "Sandusky County",
                "F": "Franklin County",
                "G": "Licking County",
                "H": "Franklin County",
                "I": "Greene County",
                "J": "Montgomery County",
            }[group],
            fine_location=group,
            enrichment="viral",
            method={
                "A": "AB",
                "B": "AB",
                "C": "C",
                "D": "D",
                "E": "EFGH",
                "F": "EFGH",
                "G": "EFGH",
                "H": "EFGH",
                "I": "IJ",
                "J": "IJ",
            }[group],
        )
    elif paper == "brinch":
        library,sample,location,date = record 
        return sample, SampleAttributes(
            date=date,
            country="Denmark",
            location="Copenhagen",
            fine_location=location,
        )
    else:
        assert False


import pprint

SampleCounts = dict[TaxID, dict[Sample, int]]

metadata_bioprojects = {} 
metadata_samples = {}
sample_counts = defaultdict(dict)
for paper, bioprojects in target_bioprojects.items():
    for bioproject in bioprojects:
        samples = []
        with open (os.path.join(BIOPROJECTS_DIR, bioproject, "sample-metadata.csv")) as inf:
            for i, record in enumerate(csv.reader(inf)):
                if i == 0:
                    continue
                sample, sample_attributes = parse_metadata(record, paper)
                samples.append(sample)
                metadata_samples[sample] = sample_attributes
        metadata_bioprojects[bioproject] = samples
        with open (os.path.join(BIOPROJECTS_DIR, bioproject, "hv_clade_counts.tsv")) as inf:
            for i, row in enumerate(inf):
                if i == 0:
                    continue
                taxid, name, rank, parent_taxid, sample, n_reads_direct, n_reads_clade = row.rstrip("\n").split("\t")
                taxid = int(taxid)
                n_reads_direct = int(n_reads_direct)
                if n_reads_direct:
                    sample_counts[taxid][sample] = n_reads_direct
        with open(os.path.join(BIOPROJECTS_DIR, bioproject, "qc_basic_stats.tsv")) as inf:
            for i, row in enumerate(inf):
                row = row.rstrip("\n").split("\t")

                if i == 0:
                    cols = row 
                    continue
                 
                metadata_samples[row[cols.index("sample")]].reads = int(row[cols.index("n_read_pairs")])


                
def load_tax_tree() -> Tree[TaxID]:
    with open("human_virus_tree-2022-12.json") as inf:
        data = json.load(inf)
    return Tree.tree_from_list(data).map(lambda x: TaxID(int(x)))


def make_count_tree(
    taxtree: Tree[TaxID], sample_counts: SampleCounts
) -> Tree[tuple[TaxID, Counter[Sample]]]:
    return taxtree.map(
        lambda taxid: (taxid, Counter(sample_counts[taxid]))
        if taxid in sample_counts
        else (taxid, Counter()),
    )


def count_reads(
    taxtree: Tree[TaxID] | None, sample_counts: SampleCounts
) -> Counter[Sample]:
    if taxtree is None:
        return Counter()
    count_tree = make_count_tree(taxtree, sample_counts)
    return sum(
        (elem.data[1] for elem in count_tree),
        start=Counter(),
    )


@dataclass
class MGSData:
    bioprojects: dict[BioProject, list[Sample]]
    sample_attrs: dict[Sample, SampleAttributes]
    read_counts: SampleCounts
    tax_tree: Tree[TaxID]

    @staticmethod
    def from_repo(
    ):
        return MGSData(
            bioprojects=metadata_bioprojects,
            sample_attrs=metadata_samples,
            read_counts=sample_counts,
            tax_tree=load_tax_tree(),
        )

    def sample_attributes(
        self, bioproject: BioProject, enrichment: Optional[Enrichment] = None
    ) -> dict[Sample, SampleAttributes]:
        samples = {
            s: self.sample_attrs[s] for s in self.bioprojects[bioproject]
        }
        if enrichment:
            return {
                s: attrs
                for s, attrs in samples.items()
                if attrs.enrichment == enrichment
            }
        else:
            return samples

    def total_reads(self, bioproject: BioProject) -> dict[Sample, int]:
        return {
            s: self.sample_attrs[s].reads for s in self.bioprojects[bioproject]
        }

    def viral_reads(
        self, bioproject: BioProject, taxids: Iterable[TaxID]
    ) -> dict[Sample, int]:
        viral_counts_by_taxid = {
            taxid: count_reads(self.tax_tree[taxid], self.read_counts)
            for taxid in taxids
        }
        return {
            s: sum(viral_counts_by_taxid[taxid][s] for taxid in taxids)
            for s in self.bioprojects[bioproject]
        }
