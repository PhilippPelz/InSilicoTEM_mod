# Modified version of InSilicoTEM v0.93

This code was used in part for generating the figures presented in [Low-dose cryo electron ptychography via non-convex Bayesian optimization](https://arxiv.org/abs/1702.05732).

It is a modified version of [InSilicoTEM v0.93](https://doi.org/10.1016/j.jsb.2013.05.008). The authors did not specify a license in their release, so I cannot include a license either. 

Changes:

-custom pdb reader
-schreibt zwischendurch die arrays, die ich fuer die ptycho-sim brauche, in hdf5 files raus
-K2 mtf und dqe

Please cite the article 

    @article{VULOVIC201319,
    title = "Image formation modeling in cryo-electron microscopy",
    journal = "Journal of Structural Biology",
    volume = "183",
    number = "1",
    pages = "19 - 32",
    year = "2013",
    note = "",
    issn = "1047-8477",
    doi = "http://dx.doi.org/10.1016/j.jsb.2013.05.008",
    url = "http://www.sciencedirect.com/science/article/pii/S1047847713001226",
    author = "Miloš Vulović and Raimond B.G. Ravelli and Lucas J. van Vliet and Abraham J. Koster and Ivan Lazić and Uwe Lücken and Hans Rullgård and Ozan Öktem and Bernd Rieger",
    }

if you use this code. If you use one of the new file starting with InSilicTEM_1ryp*, InSilicTEM_4hhb*, or InSilicTEM_4V5F*, please also consider citing 

    @article{Pelz2017,
    title = "Low-dose cryo electron ptychography via non-convex Bayesian optimization",
    journal = "Scientific Reports (accepted)",
    author = "Philipp Michael Pelz and Wen Xuan Qiu and Robert Bücker and Günther Kassier and R.J. Dwayne Miller",
    }
