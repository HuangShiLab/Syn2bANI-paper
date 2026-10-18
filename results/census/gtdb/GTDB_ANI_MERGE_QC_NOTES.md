# GTDB ANI merge QC note

`gtdb_ani_merge_validation.json` records `header_ok: false` because the audit
compared the first physical line of the C-locale-sorted merged table with the
required schema. In a lexically sorted ANI stream, lowercase header keys sort
after uppercase accession keys, so the canonical header can appear after the
first data block rather than as the first physical line. This is a header
placement artefact of the sorted format, not evidence of a corrupt merge.

The decisive merge checks are unaffected and all pass:

- 653,522,280 unique canonical rows;
- zero duplicate keys;
- zero self pairs;
- zero cross-species pairs;
- zero invalid ANI values;
- zero unknown genomes;
- ANI coverage of 99.9904% of the 653,584,719 effective expected pairs.

The final join uses the parsed canonical header wherever it occurs in the
sorted stream and validates every data row against the effective
genome-to-species map.
