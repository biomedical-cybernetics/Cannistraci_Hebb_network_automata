MATLAB code to compute the Cannistraci-Hebb (CH) network automata scores
for network links considering paths of length two (L2) or three (L3)


### REFERENCE ###

"Local-community network automata modelling based on length-three-paths
for prediction of complex network structures in protein interactomes, food webs and more"
A. Muscoloni, I. Abdelhamid, C. V. Cannistraci, bioRxiv, 2018
https://doi.org/10.1101/346916


### INPUT ###

x - adjacency matrix of the network;
    the network is considered unweighted, undirected and zero-diagonal

L - [optional] integer to indicate the path length to compute:
    0 -> compute both CH2-L2 and CH2-L3
    2 -> compute only CH2-L2
    3 -> compute only CH2-L3
    if not given or empty, the option 0 is considered

w - [optional] 2-columns matrix (id1,id2) indicating the links for which the score should be calculated;
    if not given or empty, the scores for all the missing links are computed

par - [optional] 1 or 0 to indicate whether the function should use parallel computation or not;
    if not given or empty, parallel computation is used


### OUTPUT ###

scores - 3-columns or 4-columns matrix depending on the input parameter L:
    L=0 -> 4-columns matrix containing the values (id1,id2,score_CH2_L2,score_CH2_L3)
    L=2 -> 3-columns matrix containing the values (id1,id2,score_CH2_L2)
    L=3 -> 3-columns matrix containing the values (id1,id2,score_CH2_L3)
 

### CONTACT ###

For any problem, please contact:
Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com