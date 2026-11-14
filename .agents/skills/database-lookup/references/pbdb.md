# Paleobiology Database (PBDB) API v1.2 Reference

## Overview
Community-curated authority for fossil occurrences, taxonomic opinions, stratigraphic ranges (First and Last Appearance Data, FAD/LAD in Ma), paleoecological traits, and geological time intervals across the Phanerozoic and Precambrian.

## Base URL
```
https://paleobiodb.org/data1.2/
```

## Authentication
- **No API key required** for read endpoints (`taxa/*`, `occs/*`, `intervals/*`, `strata/*`, `specs/*`, `refs/*`).
- Identify automated clients via `User-Agent` header and respect fair-use rate limits (≤ 3 req/sec; use single-call filters rather than per-species loops when possible).

## Key Endpoints

### 1. Single Taxon & Stratigraphic Range (`taxa/single.json`)
Retrieves taxonomic rank, status, extant/extinct flag (`is_extant`), and stratigraphic first/last appearance bounds (`firstapp_max_ma`, `firstapp_min_ma`, `lastapp_max_ma`, `lastapp_min_ma` when `show=app`).

```
GET https://paleobiodb.org/data1.2/taxa/single.json?name=TAXON&show=attr,app,parent,ecotaph
```

| Param | Description |
|-------|-------------|
| `name` | Scientific name (e.g., `Tyrannosaurus rex`, `Archosauria`, `Spinosaurus`) |
| `id` | PBDB taxon ID (e.g., `txn:38613` or integer `38613`) |
| `show` | `app` (FAD/LAD appearance ages in Ma), `attr` (author/year), `parent`, `ecotaph` (environment, diet, motility, life habit), `size` (extant/total subtaxa count) |
| `vocab` | `pbdb` (human-readable field names: `max_ma`, `min_ma`, `taxon_name`) or `com` (compact default: `fea`, `lla`, `nam`) |

**Example (compact vs `vocab=pbdb`):**
```bash
curl --fail-with-body --silent --show-error \
  "https://paleobiodb.org/data1.2/taxa/single.json?name=Tyrannosaurus+rex&vocab=pbdb&show=app,ecotaph"
```
Key fields in `records[0]` (`vocab=pbdb`):
- `orig_no` / `taxon_no`: PBDB identifier
- `taxon_name`, `taxon_rank`, `is_extant` (`"extant"` or `"extinct"`)
- `firstapp_max_ma`, `firstapp_min_ma`: FAD uncertainty interval (Ma)
- `lastapp_max_ma`, `lastapp_min_ma`: LAD uncertainty interval (Ma)
- `environment`, `diet`, `motility`, `vision`

### 2. Clade Subtree & Children (`taxa/list.json`)
```
GET https://paleobiodb.org/data1.2/taxa/list.json?base_name=Archosauria&rank=family&vocab=pbdb&show=app
```

| Param | Description |
|-------|-------------|
| `base_name` | Clade name to traverse downward |
| `rel` | `children` (direct), `all_children` (subtree), `synonyms`, `common_ancestor` |
| `rank` | Filter by rank (`species`, `genus`, `family`, `order`, `class`) |
| `extant` | `true` or `false` |

### 3. Fossil Occurrences (`occs/list.json`)
Retrieves geo-referenced and stratigraphically dated fossil occurrences for PyRate, Fossilized Birth-Death (FBD) sampling rate ($\psi$), or biogeographic reconstruction (DEC/BayArea).

```
GET https://paleobiodb.org/data1.2/occs/list.json?base_name=Spinosauridae&vocab=pbdb&show=coords,phylo,strat,time&limit=500
```

| Param | Description |
|-------|-------------|
| `base_name` | Taxon or clade name (includes all synonyms and subtaxa automatically) |
| `interval` | Geological interval name (e.g., `Cretaceous`, `Cenomanian`, `Triassic`) |
| `max_ma`, `min_ma` | Numerical age bounds in Ma |
| `cc` | ISO2 country code or continent code (`AF`, `SA`, `NA`, `EU`, `AS`, `OC`, `AN`) |
| `show` | `coords` (modern `lng`/`lat`), `paleoloc` (paleocoordinates `paleolng`/`paleolat`), `phylo` (taxonomy), `strat` (formation, member), `time` (`max_ma`, `min_ma`, `early_interval`, `late_interval`) |
| `limit` | Max records (`all` or integer; use `limit=1` with `count=true` first to check total) |

### 4. Geological Time Intervals (`intervals/list.json`)
```
GET https://paleobiodb.org/data1.2/intervals/list.json?scale=1&vocab=pbdb
```
Returns ICS International Chronostratigraphic Chart intervals (`interval_name`, `early_age`, `late_age` in Ma).

## Filter & Completeness Hazards
1. **Default vocabulary (`com`) uses 3-letter keys** (`nam`, `fea`, `lla`, `eag`, `lag`). Always pass `vocab=pbdb` for self-documenting keys (`taxon_name`, `firstapp_max_ma`, `lastapp_min_ma`, `max_ma`, `min_ma`).
2. **Occurrence age is an interval, not a point**: Every fossil occurrence has `[min_ma, max_ma]` (the boundaries of the stratigraphic bin). Never average `(max_ma + min_ma) / 2` as an exact tip date in BEAST2/RevBayes; use a `uniform(min_ma, max_ma)` tip-age prior in FBD or PyRate's randomized age-bin resampling.
3. **Node calibration vs. First Appearance**: The oldest known fossil of a clade (`firstapp_max_ma` / `firstapp_min_ma`) provides only a **hard minimum bound** (`min_ma` of that specimen's stratigraphic unit) for the crown/total-group divergence, never a hard maximum.
