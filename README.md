# Open data of DNS tunneling

The repository contains datasets for DNS tunneling detection. Reference paper is [1]. 

The dataset is composed by statistics collected on DNS messages (without and with tunnel): average, variance, skewness and kurtosis of interarrival times [s], query size [B], answer size [B]. The last column of the dataset contains a boolean corresponding to the presence of tunnel '1' or not '0'. 

The applications tunneled in regular DNS traffic are SSH and p2p through DNS2TCP tool. Regular traffic has been collected @ CNR premises, IEIIT, Genova, Italy. The first half of each dataset is with label 0 and the second with 1. The tunnel traffic is composed of 90% of legitimate traffic and 10% of tunneled application. (No tunnel means 100% of legitimate DNS traffic). Details of DNS network, monitoring and features extraction are reported in [2].

The rule-based models for detection are also reported as well as the Matlab code for data visualization. Rules are derived through the logic learning machine in www.rulex.ai. Other approaches, such as decision trees and rule-extraction from deep neural networks, are studied in the reference paper.

As to main.cpp, the following reccomandations hold. It is the code for dataset generation. It refers to dns2tcpssh.txt and dns2tcpp2p.txt where the two applications are tunneled inside regular DNS. The code creates a mix of DNS (label 0 in the datasets) with tunnels (for ssh and p2p separately, label 1 in the datasets). Just replace dns2tcpp2p.txt with dns2tcpssh.txt in the main code to swhitch to the desired application. The "contatoreMix" variable in the code should be set as follows in order to obtain different mixtures of hidden application inside regular DNS: 1 corresponding to 50% - 9 10% - 99 1% - 999 0.1% - 5 20%.

For larger datasets or any issue: maurizio.mongelli@cnr.it.

[1] G. De Bernardi, G. Battista Gaggero, F. Patrone, M. Mongelli, M. Marchese, "Rule-based eXplainable Autoencoder for DNS Tunneling Detection," Computer Communication, submitted.

[2] Aiello, M., Mongelli, M., and Papaleo, G. (2015) DNS tunneling detection through statistical fingerprints of protocol messages and machine learning. Int. J. Commun. Syst., 28: 1987–2002. doi: 10.1002/dac.2836.
