# ARU data and N-mixture models

ARU recordings are often divided into standardized clips, but those clips are
not automatically valid N-mixture "visits."

For a standard N-mixture model, each site needs repeated counts while the local
population is closed, and counts must be conditionally independent given the
latent abundance $N_i$ and detection probability $p_{ij}$. A defensible ARU
design is therefore usually several fixed-duration recordings on different days
(or sufficiently separated survey occasions) within a short breeding-season
window. Time of day, date, weather, recorder, and clip length can be modeled in
$p$.

Splitting one continuous 30- or 60-minute recording into 3- or 5-minute pieces
can be useful, but it is generally better interpreted as within-survey data for
time-removal or availability models. The same vocal individual drives strong
serial dependence across adjacent clips, so treating them as independent
N-mixture visits can badly overstate information and make abundance estimates
fragile. With a single ARU, a "count" is also usually a count of vocalizations
or detection events, not reliably distinct individuals; that makes the usual
binomial observation model especially questionable.

For ARU data, use time-removal within a recording, and reserve N-mixture for
genuinely repeated, standardized recordings across days.

## References

- Furnas, B. J., and R. L. Callas (2015). Using automated recorders and occupancy models to monitor common forest birds across a large geographic region. *Journal of Wildlife Management*, 79, 325-337. [DOI](https://doi.org/10.1002/jwmg.821)

Furnas and Callas provide an ARU design example: 5-minute recordings at three morning times, repeated for three consecutive days. They use single-season occupancy, not N-mixture.

- Thompson, S. J., C. M. Handel, and L. B. McNew (2017). Autonomous acoustic recorders reveal complex patterns in avian detection probability. *Journal of Wildlife Management*, 81, 1228-1241. [DOI](https://doi.org/10.1002/jwmg.21285)

Thompson et al. use time-removal methods on 10-minute ARU recordings to estimate availability and discuss short recording intervals.

## Recent ARU Literature (2022-2024)

This review is independent of the project material above. Recent applied ARU
studies more often use repeated-visit occupancy models than conventional
N-mixture models. They convert recordings to detection histories, generally
require at least three repeat samples per site, and account for recorder,
recording, and classification error in the observation model.

The most defensible repeats are standardized recordings from separate survey
occasions, often on multiple days during a short closure period. ARUs make this
design practical because they can record the same site without repeated field
visits. Some studies also divide recordings into fixed-length segments, but
these segments should be regarded as repeated detection opportunities only after
considering temporal dependence. Consecutive clips dominated by the same vocal
individual are not automatically independent visits.

The distinction matters for the estimand. An occupancy model describes the
probability that a species occupies a defined site or sampling area; it does not
require assigning calls to individuals. An N-mixture model requires repeated
counts that can plausibly be treated as binomial detections of a closed local
abundance. That requirement is harder to justify for ARU call or event counts,
so occupancy is currently the more common formulation.

### Recent examples

- Cole, J. S., N. L. Michel, S. A. Emerson, and R. B. Siegel (2022). *Automated bird sound classifications of long-duration recordings produce occupancy model outputs similar to manually annotated data.* *Ornithological Applications*, 124(2), duac003. [DOI](https://doi.org/10.1093/ornithapp/duac003)

The clearest bird example. The study states that occupancy estimation needs a minimum of three repeat samples, analyzes 9-minute recording segments as repeated samples, and incorporates validated false-positive and true-positive rates in the occupancy model. It also compares 9- and 87-minute automated annotations.

- Wood, C. M., and M. Z. Peery (2022). *What does 'occupancy' mean in passive acoustic surveys?* *Ibis*, 164(4), 1295-1300. [DOI](https://doi.org/10.1111/ibi.13092)

This interpretation paper distinguishes recorders placed at known-use areas from recorders on a random grid, and argues that the biological meaning of occupancy, detection, colonization, and extinction depends on that design choice.

- Wood, C. M., A. Barceinas Cruz, and S. Kahl (2023). *Pairing a user-friendly machine-learning animal sound detector with passive acoustic surveys for occupancy modeling of an endangered primate.* *American Journal of Primatology*, 85(8), e23507. [DOI](https://doi.org/10.1002/ajp.23507)

Demonstrates a single-season occupancy model based on passive acoustic detections of Yucatan black howler monkeys. It is useful beyond birds because it shows the now-common workflow of classifier output followed by an explicit detection model.

- James, R., J. R. Bennett, S. Wilson, et al. (2024). *Modelling the occupancy of two bird species of conservation concern in a managed Acadian Forest landscape: Applications for forest management.* *Forest Ecology and Management*, 555, 121725. [DOI](https://doi.org/10.1016/j.foreco.2024.121725)

A direct example of repeated ARU monitoring: recordings were collected across multiple days during the breeding season and used in occupancy models for Canada Warbler and Olive-sided Flycatcher habitat associations.
