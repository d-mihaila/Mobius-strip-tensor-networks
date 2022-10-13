# Mobius-strip-tensor-networks

This code is computing the Hamiltonian of a superconducting Mobius strip. Using a DMRG algorithm via the pytenet library we calculate the minimum energies for variouse circuit sizes. 
Our goal is to analize the coherence properties of the 0-π qubit as it aproaces near degenerate ground state. 

This research is based on the article 'Spectrum and coherence properties of the current-mirror qubit' and we used the quantum tensor network library by Professor Christian Mendl cmendl/pytenet as referenced below.

Note: Our resulting graph does not have the expeted shape and this, we believe, is due to the periodic boundary conditions of our system, for which the Pytenet library will be adjusted. 


# References 
1. Weiss, D. K., Andy CY Li, D. G. Ferguson, and Jens Koch. "Spectrum and coherence properties of the current-mirror qubit." Physical Review B 100, no. 22 (2019): 224507.
   
2. Mendl, Christian B. "PyTeNet: A concise Python implementation of quantum tensor network algorithms." Journal of Open Source Software 3, no. 30 (2018): 948.
