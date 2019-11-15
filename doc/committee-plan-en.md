## Supervisory committee plan - Paragraph by paragraph

### Introduction

1. Local contribution to beta-diversity (LCBD) calculation methods allow for the
   identification of sites with exceptional biodiversity.
   However, such methods are currently restricted to discontinuous spatial scales.

2. Species distribution models (SDM) already allow us to make predictions on continuous
   spatial scales.
   Hence, they could possibly be used for LCBD calculations.

3. The BIOCLIM and eBird databases are frequently increasingly used in SDMs, hence they could
   also be used to calculate LCBDs on continuous spatial scales.

4. This approach would be especially useful in regions with only sparse sampling.

5. The method could also be applied on IPCC climate change scenarios in order to identify
   sites where beta-diversity changes will be most important.

### Methods

1. Data Sampling:
   We collected the data available in the databases for North America and for the Warblers
   family.

2. Data Manipulation : We converted raw observation data into presence-absence per pixel at a
   resolution of 10 arc-minutes.
   (We also applied Hellinger transformation, although the most appropriate method remains to
   be determined)

3. SDM : We used the BIOCLIM method to predict species distributions based on the quantiles
   of environmental variables at observed sites.

4. LCBD : We calculated the LCBD statistics through the total variance of the matrix Y,
   instead of though a dissimilarity matrix.

5. Prediction validity:
   The exact way of testing the validity of the predictions remains to be determined.

6. Alternative methods : Other methods could possibly outperform BIOCLIM for the predictions,
   including MAXENT, Random Forests, neural networks, Joint SDMs, etc.

7. Climate change scenarios & temporal beta diversity : We aim to apply our method to
   environtal conditions from climate change scenarios, and then use temporal beta-diversity
   indices (TBI) to identify the sites with the most important changes.

8. Preliminary results : Our preliminary results mainly compare raw data statistics to
   prediction statistics.
   The figures will focus on :
   - Individual species distribution
   - Richness
   - LCBD
   - Richness-LCBD relationship
