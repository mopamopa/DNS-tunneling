# Open data of DNS tunneling

The repository contains datasets for DNS tunneling detection. 

Reference paper is "Rule-based eXplainable Autoencoder for DNS Tunneling Detection," Computer Communication journal, Elsevier. By G. De Bernardi, G. Battista Gaggero, F. Patrone, M. Mongelli, M. Marchese.

The two applications tunneled in regular DNS traffic are SSH and p2p through DNS2TCP tool. Regular traffic has been collected @ CNR premises, IEIIT, Genova, Italy. 

The dataset is composed by statistics collected on DNS messages (without and with tunnel): average, variance, skewness and kurtosis of interarrival times [s], query size [B], answer size [B]. The last column of the dataset contains a boolean corresponding to the presence of tunnel '1' or not '0'. The first half is with label 0 and the second with 1.

The rule-based models for detection are also reported as well as the Matlab code for data visualization. Rules are derived through the logic learning machine in www.rulex.ai. Other approaches, such as decision trees and rule-extraction from deep neural networks are studied in the paper.

maurizio.mongelli@cnr.it.
