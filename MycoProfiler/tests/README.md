# Test data

Ten synthetic read pairs (`sample_R1.fastq` / `sample_R2.fastq`) and the same ten
reads unpaired (`sample_single.fastq`). Read names carry a tag the mock classifier
keys on:

| tag        | pairs | represents                                                |
|------------|-------|-----------------------------------------------------------|
| `BACTERIA` |     5 | classified at stage 1, removed                            |
| `FUNGUS`   |     3 | unclassified at stage 1, then classified at stage 2       |
| `NOVEL`    |     2 | unclassified at both stages                               |

These are not real sequences and carry no biological meaning; they exist so the
test suite can assert exact counts at every stage without a multi-gigabyte database.
