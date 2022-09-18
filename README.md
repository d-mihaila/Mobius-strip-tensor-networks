# Mobius-strip-tensor-networks

This code is computing the Hamiltonian of a superconducting Mobius strip. Using a DMRG algorithm via the pytenet library we calculate the minimum energies for variouse circuit sizes. 
Our goal is to analize the coherence properties of the 0-π qubit as it aproaces near degenerate ground state. 

This research is based on the article 'Spectrum and coherence properties of the current-mirror qubit' and we used the quantum tensor network library by Professor Christian Mendl cmendl/pytenet as referenced below.

Note: Our resulting gaph does not have the expeted shape and this, we believe, is due to the periodic boundary conditions of our system, for which the Pytenet library will be adjusted. 


# References 
1. D. K. Weiss, Andy C. Y. Li, D. G. Ferguson, Jens Koch
   | Spectrum and coherence properties of the current-mirror qubit
   | arXiv:1908.04615
   
2. @ARTICLE{pytenet,
  author = {Mendl, C. B.},
  title = {PyTeNet: A concise Python implementation of quantum tensor network algorithms},
  journal = {Journal of Open Source Software},
  year = {2018},
  volume = {3},
  number = {30},
  pages = {948},
  doi = {10.21105/joss.00948},
}
