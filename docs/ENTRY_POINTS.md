# Study entry points

Start with `make demo` when you do not have study data. The entries below are a
selected navigation index, not a declaration that every entry runs with public
inputs alone. See [data requirements](REPRODUCIBILITY.md).

| Topic | Register ID | State / compute | Input availability |
|---|---|---|---|
| Build/check the multi-sensor cube | `cube_build` | current / mac | Public sources plus restricted field inputs |
| Download RTC stacks | `rtc_download` | current / mac | Provider access and configured storage |
| Inspect Sentinel-1 coverage | `s1_extension_candidates` | current / mac | Public catalogue/provider access |
| Coregister the SLC pilot | `slc_cslc_x101` | finished / colab | Earthdata, cloud SAR environment and study geometry |
| Score the baseline chain | `baseline_x099` | current / mac | Private cube, truth, parameters and pre-registration |
| Compare short pairs with the laser | `field_p6_short_pairs_x065` | finished / mac | Private truth and radar inputs |
| Evaluate winter model selection | `nested_winter_x135` | finished / mac | Private frozen study artifacts and truth |
| Inspect laser datum/fill sensitivity | `laser_qc_x147` | current / mac | Private laser record and study artifacts |
| Analyze spatial timing | `p6_motion_footprint_x133` | current / mac | Private study inputs |
| Build the short-pair comparison figure | `article_fig_f4` | open / mac | Private result tables |
| Export an atlas data package | `build_web_atlas` | current / mac | Private study products and output workspace |

The complete [register](../config/register.yaml) is authoritative. Query it with:

```bash
make list
make list KIND=pipeline
make list TICKET=X-065
```

`current` means a maintained entry, `finished` means a recorded study execution,
and `open` means planned work. An archived entry has a non-HEAD `code_ref`; inspect
that ref rather than expecting the script in the current tree. `mac` denotes the
local/light-compute route; `colab` entries require the documented cloud setup.

Ticket and try IDs are stable cross-references to the study record. The private
ledger is not required to read the public implementation. The public register
contains technical inputs/outputs and questions; it does not include private
field tables or the private decision log.

Run study entries through `make run E=<id>` after configuring their inputs.
See the [figure/table map](MANUSCRIPT_MAP.md) for manuscript navigation.
