# Codebook — `wtp_data.csv`

One row per respondent, N = 71. Rows are in random order (see *De-identification*).

| Column | Type | Unit | Meaning |
|---|---|---|---|
| `priv_price` | integer | TSh | Willingness to pay for tree planting on **private** land (`binafsi`). The highest bid the respondent accepted on the bid ladder, inclusive of uncertain (`?`) answers. |
| `pub_price` | integer | TSh | The same quantity for **community** land (`jamii`). |
| `sex` | character | — | `M` / `F`. |
| `ms` | character | — | Marital status. |
| `age` | integer | years | Age in 2023, derived as `2023 - birth_rev` (birth year only; no date of birth was recorded). |
| `ed_rev` | integer | years | Completed years of education, revised. |
| `wealth` | numeric | TSh | Total household wealth in 2023, corrected (see *Provenance*). |

The models do not use these columns directly. `code/prepare_data.R` builds the
Stan inputs from them: `priv` and `pub` are `log1p` of the two prices, `wealth`
is standardized `log1p`, `age` is standardized, and `sex` and `ed_rev` become
1-based integer indices.

## Elicitation

Willingness to pay was elicited with a bid ladder in steps of 1,000 TSh. The stored
value is the *high* estimate — the highest accepted bid counting uncertain
answers as acceptances. The dataset's original author notes that the high and
low variants were compared and "results looked the same".

## Provenance and corrections

Transcribed from the `info_updated` sheet of the source workbook, which is the
only place this documentation previously existed.

**1 April 2025.** The revised dataset replaced `totalwealth2018` with
`totalwealth2023`; added individuals whose IDs had been uncertain in the first
dataset; and replaced missing demographic data. Typos were corrected in the
community-land low-bid column, where six respondents had `1000` entered instead
of `100`.

**29 April 2025.** `totalwealth2023` was replaced with a corrected variant. For
polygynously married men, wealth had been attributed on the basis of the first
sharing unit in which they were recorded, rather than the average across both
of their sharing units; this changed the values for two respondents. An error in
a third respondent's wealth was also corrected. `wealth` in this release is that
corrected variable.

## De-identification

This file is derived from a source workbook that contains respondent names and
panel identifiers. Removed before release:

- `name` and `named` — respondent full names.
- `PE` — the individual respondent identifier, which links to the ENDOW panel
  study and whose prefix encodes the sharing unit, so that two respondents from
  one household could be recognised as such.
- `PESU2018` — the sharing-unit (household) identifier.
- One respondent with no revised birth year was dropped, following the analysis
  script, taking 72 rows to 71.

Row order was randomised under a fixed seed. This matters: in the source
workbook rows are sorted by `PE`, and `PE` order is alphabetical by respondent
name, so row position alone carried identifying information. The models are
exchangeable across rows, so the permutation does not affect any result.

The released file is **pseudonymous, not anonymous.** See the *Data* section of
the top-level README for the residual re-identification risk.
