# Epidemiological interventions illustrations engine

Kucharavy Andrei 2020 all rights reserved - MIT license.
Generated illustrations are CC BY-SA.

## Disclaimer
This engine is written for illustrative purposes only.

The assumptions of the simulations are overly simplifying and it is not meant to be used for any sort of predictions - just illustrate some basic concepts. This engine cannot be used to simulate actual outbreaks without major modifications.

You will need a basic scientific python 3.8 installation to run this. Anaconda one is recommended.


## Assumptions
We use the following assumptions about the disease:
 - 5 days, non-infective incubation period
 - 48 hours pre-symptomatic infectious period
 - 14 days symptomatic infections period
 - after symptoms subside, the population acquires immunity
 - All of the above can be adjusted using the `state_duration` list

We use the following assumptions about the transmission:
 - any points within a close circle ("sure contaminations") will be contaminated as soon as an
 individual becomes infectious. It can be close coworkers/household members/hobby activity contacts
 - at most 1 individual in a larger circle ("possible contaminations") can be contaminated per
 day. The chance of contamination is proportional to the distance to the contagious individual.
 There is a chance that on a given day no contamination occurs that is set as a parameter.
 - All of the above can be adjusted using the `limit_distances` list

We use the following assumptions about the measures:
 - The quarantine is not perfect, but reduces the radius of "sure contaminations" and "possible
 contaminations" circles as well as raises the probability no contamination occurs on a given day in
 "possible contaminations" circle (eg transmissions within the household, during store visits, ...)
 - The degree of efficiency of quarantine is controlled by the `quarantine_power` parameter
 - The adherence to all the measures (quarantine, self-isolation, contact tracing & isolation...) is
 total
 - Contact tracing can vary between two states: "from memory only", where only the "sure
 contaminations" are isolated and "total tracking", where in addition all the "possible
 contaminations" are isolated as well.
 - The degree of efficiency of the contact tracing power is controlled by the
 `contact_tracing_power` parameter.

## R_t calculation
R_t is calculated as the average number of individuals infected by a single infectious
individual across all the individuals that are still infections

## Sample runs: 
| scenario | relevant parameters | total infections | symptomatic peak | days until end |
|---|---|---|---|---|
|no intervention (herd immunity) | default  |  896 |  282 |  107 |
|hardcore quarantine till the end |  `quarantine_day = 20`, `quarantine_power = 4`  | 97 |  80 |  63 |
|too soft quarantine   | `quarantine_day = 20`, `quarantine_power = 2`  | 846 | 129  | 211  |
|too soft quarantine for 6 weeks   | `quarantine_day = 20`,  `quarantine_end_day  = 40`, `quarantine_power = 4`  |  894 |  288 | 132  |
|self-isolation alone   | `self_isolation_day = 20`, `quarantine_power = 4`  | 896  | 285  | 106  |
|hardcore contact tracing   | `contact_tracing_day = 20`, `quarantine_power = 4`, `contact_tracing_power = 1.` | 167  | 80  | 77  |
|minimal contact tracing   | `contact_tracing_day = 20`, `quarantine_power = 4`, `contact_tracing_power = 0.` | 896  | 264  |  122 |
|90% vaccination_rate   | `fraction_initially_immune = 0.9` | 3  | 3  |  28 |
|70% vaccination_rate   | `fraction_initially_immune = 0.7` | 301  | 81  |  122 |
