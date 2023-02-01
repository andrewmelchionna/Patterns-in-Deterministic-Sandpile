# Patterns-in-Deterministic-Sandpile
The Abelian Sandpile Model is a system of interacting particles. We explore the algebraic structure of the model and the patterns that it generates.

# Deterministic Sandpile on Ellipse
In the Abelian Sandpile Model, particles live on the vertices of a graph, and the particle configuration updates in discrete timesteps according to the following rule.
If at any given time, a vertex is hosting at least as many particles as its degree, it "topples" and sends one particle to each of its neighbors. 
If we designate a vertex in our graph as the "sink" vertex (which never topples; particles falling onto the sink vertex stay there forever), this system will eventually stabilize.

This project takes elliptical subsets of the square grid as the underlying graph, with the boundary of the ellipse acting as the sink vertices. By running the sandpile model on 
ellipses with a certain shape, we see the sandpile "resonate" and produce periodic patterns. These experiments demonstrate the interesting fractal geometry contained implicitly
in the model's simple rules: https://arxiv.org/abs/1208.4839 !
