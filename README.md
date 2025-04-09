# CSC 448 - Spring 2025 #
## [Lab 0](https://github.com/Pilon2020/CSC448/tree/main/Lab0)
Lab 0 was just a lab to get everything set up. Importing key libraries into workspace for future labs.  
- Imported Numpy, Scipy and sklean into workspace

### [Smith Waterman Algorithm](https://github.com/Pilon2020/CSC448/tree/main/Lab0/SmithWaterman.py)
Working on getting a my own version of the Smith-Waterman algorithm working. Much of the information which I am checking against is coming from the [Wikipedia Article](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) on the topic. I found some previous code I found on [github](https://github.com/slavianap/Smith-Waterman-Algorithm/blob/master/Script.py) had published that I used as an aid whenever I got stuck.  
  
  
I have been able to get it to work with a basic gap penalty but as soon as I add in an Affine penalty in it breaks everything.

## [Project 1](https://github.com/Pilon2020/CSC448/tree/main/Lab1)  
Project 1 - Sequences and Evolutionary Trees  
Using the Smith-Waterman Algorithm I worked on last time, while adding the BLOSUM62 weighting system in allows me to accurately score the similarities between protien strings.  
- Scores are normalized by dividing the compared Smith-Waterman (SW) score by the largest SW score produced when comparing the same string: SW(S<sub>i</sub>,S<sub>j</sub>) / max(SW(S<sub>i</sub>,S<sub>i</sub>),SW(S<sub>j</sub>,S<sub>j</sub>))
- The normalized scores were then stored in an array.  
Still need to implament a clustering algorithm to construct a phylogenic tree.
