### About
The code was written as part of my master thesis project, where I used neural networks to obtain optimal hedging strategies for a variety of option types in the Black-Scholes and Heston model. The network is implemented in Jupyter Notebook.

The network structure is based on Mikko Pakkanens code and slides, and the Deep Hedging article by Buehler et al. https://arxiv.org/abs/1802.03042. 
Data from the Heston model was generated using the Broadie-Kaya exact simulation scheme, implemented in C++. 
The C++ code is based on Dave Roberts code found at: https://github.com/daleroberts/heston-qmc. 
It requires the package Boost found at https://www.boost.org/

By Nicolai Munch Kofoed. 

Mail: Nicolaisteenfos@gmail.com
