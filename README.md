# Open data of DNS tunneling

The repository contains datasets for DNS tunneling detection as to paper Rule-based eXplainable Autoencoder for DNS Tunneling Detection, Computer Communication journal, Elsevier. By G. De Bernardi, G. Battista Gaggero, F. Patrone, M. Mongellia, M. Marchese.
The two applications tunneled are SSH and p2p. The rule-based models for detection are also reported.
The dataset is composed by statistics collected on DNS messages (without and with tunnel): average, variance, skewness and kurtosis of interarrival times [s], query size [B], answer size [B]. The last column of the dataset contains a boolean corresponding to the presence of tunnel '1' or not '0'. The first half is with label 0 and the second with 1.
maurizio.mongelli@cnr.it.
