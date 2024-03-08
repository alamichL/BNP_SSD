# BNP-SSD
This repository presents the code used in the paper: **Species Sensitivity Distribution revisited: a Bayesian nonparametric approach**, *Louise Alamichel, Julyan Arbel, Guillaume Kon Kam King
and Igor Prünster*. 

# BNP-SSD Shiny application
A Shiny application named BNP-SSD is available at XXX. 
This application, inspired by the `shinyssdtools` Shiny application, is based on the functions of the package `BNPdensity`. 
The BNP model described in the paper is adjusted on datasets either censored or not. The fitted density is plotted as well as some goodness of fit plots. Finally, the induced optimal clustering is computed and plotted.

## Abstract
We present a novel approach to ecological risk assessment by reexamining the Species Sensitivity Distribution (SSD) method within a Bayesian nonparametric (BNP) framework. Widely mandated by environmental regulatory bodies globally, SSD has faced criticism due to its historical reliance on parametric assumptions when modeling species variability. By adopting nonparametric mixture models, we address this limitation, establishing a more statistically robust foundation for SSD. 

Our BNP approach offers several advantages, including its efficacy in handling small datasets or censored data, which are common in ecological risk assessment, and its ability to provide principled uncertainty quantification alongside simultaneous density estimation and clustering. We utilize a specific nonparametric prior as the mixing measure, chosen for its robust clustering properties—a crucial consideration given the lack of strong prior beliefs about the number of components. 

Through systematic simulation studies and analysis of real datasets, we demonstrate the superiority of our BNP-SSD over classical normal SSD and kernel density estimate SSD methods. We provide a Shiny application, BNP-SSD, making our methods available to the ecotoxicology community. 

Moreover, we exploit the inherent clustering structure of the mixture model to explore patterns in species sensitivity. Our findings underscore the effectiveness of our approach in improving ecological risk assessment methodologies.
