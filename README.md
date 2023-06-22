  <!-- badges: start -->
  [![R-CMD-check](https://github.com/ericward-noaa/bycatch/workflows/R-CMD-check/badge.svg)](https://github.com/ericward-noaa/bycatch/actions)
  <!-- badges: end -->
  
# bycatch
This is a repository for applications of Bayesian bycatch models, using Stan for estimation. These models fit Bayesian GLMs (with or without time varying parameters, and with or without covariates) to estimate bycatch rates. The package includes functions to estimate fleet-level estimates. The primary applications of these models (linked to below) are protected species bycatch (e.g. ESA listed birds, marine mammals) in the west coast groundfish fisheries off the west coast of the USA. Data compilation and analysis is done at the Northwest Fisheries Science Center (NWFSC) in Seattle. The pkgdown version of this site can be found here: [https://ericward-noaa.github.io/bycatch/](https://ericward-noaa.github.io/bycatch/)

## Citation:

Paper: 
Jannot, J.E., Ward, E.J., Somers, K.A., Feist, B.E., Good, T.P., Lawson, D., and J.V. Carretta. 2021. Using Bayesian time-series models to estimate humpback whale entanglements in the U.S. west coast sablefish pot fishery. Frontiers in Marine Science. 

Package: 
Ward, E.J. and J. Jannot. bycatch: Using Bayesian generalized linear models for estimating bycatch rates and generating fleet-level expansions.  [![DOI](https://zenodo.org/badge/85732013.svg)](https://zenodo.org/badge/latestdoi/85732013)

**See also**:  

Stohs, S.M. 2008. Predicting effort and protected species bycatch under an effort limit or take caps. American Agricultural Economics Association Annual Meeting, Orlando, FL.

Gardner, B., Sullivan, P.J., Epperly, S., and Morreale, S.J. 2008. Hierarchical modeling of bycatch rates of sea turtles in the western North Atlantic. Endangered Species Res. 5(December): 279–289. https://doi:10.3354/esr00105

Martin, S.L., S.M. Stohs, and J.E. Moore. 2015. Bayesian inference and assessment for rare‐event bycatch in marine fisheries: a drift gillnet fishery case study. Ecological Applications, 25(2):416-429. https://doi.org/10.1890/14-0059.1 

Good, T.P., J.E. Jannot, K.A. Somers, and E.J. Ward. 2022. Using Bayesian time series models to estimate bycatch of an endangered albatross. Fisheries Research,
256:106492. https://doi.org/10.1016/j.fishres.2022.106492

## Reports:

Jannot, J. E., A. Wuest, T. P. Good, K. A. Somers, V. J. Tuttle, K. E. Richerson, R. S. Shama, and J. T. McVeigh. 2021. Seabird Bycatch in U.S. West Coast Fisheries, 2002–18. U.S. Department of Commerce, NOAA Technical Memorandum NMFS-NWFSC-165. [doi: https://doi.org/10.25923/78vk-v149](https://doi.org/10.25923/78vk-v149)

Jannot, J. E., K. A. Somers, V. Tuttle, J. McVeigh, J. V. Carretta, and V. Helker. 2018. Observed and Estimated Marine Mammal Bycatch in U.S. West Coast Groundfish Fisheries, 2002-16. U.S. Department of Commerce, NWFSC Processed Report 2018-03. [doi: https://doi.org/10.25923/fkf8-0x49](https://doi.org/10.25923/fkf8-0x49)

Good, T. P., E. Ward, J. Jannot, R. Shama, N. Riley, and J. McVeigh. 2019. Observed and Estimated Bycatch of Short-tailed Albatross in U.S. West Coast Groundfish Fisheries 2002-2017. National Marine Fisheries Service, NWFSC, 2725 Montlake Blvd E., Seattle, WA 98112 [link](https://www.pcouncil.org/documents/2019/06/agenda-item-i-4-a-nmfs-report-6-observed-and-estimated-bycatch-of-short-tailed-albatross-in-u-s-west-coast-groundfish-fisheries-2016-2017-electronic-only.pdf/)

Hanson, M.B., T.P. Good, J.E. Jannot, and J. McVeigh. 2019. Estimated humpback whale bycatch in the U.S. West Coast Groundfish Fisheries 2002-2017. National Marine Fisheries Service, NWFSC, 2725 Montlake Blvd E., Seattle, WA 98112

Good, T. P., E. Ward, J. Jannot, R. Shama, N. Riley, and J. McVeigh. 2017. Observed and Estimated Bycatch of Short-tailed Albatross in U.S. West Coast Groundfish Fisheries 2014-2015. National Marine Fisheries Service, NWFSC, 2725 Montlake Blvd E., Seattle, WA 98112 [link](https://www.pcouncil.org/documents/2017/04/agenda-item-f-5-a-nmfs-report-6.pdf/)

## NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [NOAA Fisheries](https://www.fisheries.noaa.gov/)
