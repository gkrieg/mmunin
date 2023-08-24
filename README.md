# mmunin

## About the shortest hyperpath problem in systems biology

Cell signaling pathways are cornerstones of molecular and cellular biology. They underly cellular communication, govern environmental response, and their perturbation has been implicated in the cause of many diseases. While signaling pathways have classically been modeled as ordinary graphs, using directed or undirected edges to link pairs of interacting molecules, it has been shown that ordinary graphs cannot adequately represent cellular activity that involves the assembly and disassembly of protein complexes, or multiway reactions among such complexes. We instead model these pathways using directed hypergraphs, a generalization of ordinary graphs where an edge, now called a hyperedge, is directed from one set of vertices, called its tail, to another set of vertices, called its head.

Biologically, a typical cell-signaling pathway consists of membrane-bound receptors that bind to extracellular ligands, triggering intracellular cascades of reactions, culminating in the activation of transcriptional regulators and factors. Computationally, treating receptors as sources, and transcription factors as targets, finding the most efficient way to synthesize a particular transcription factor from a set of receptors maps to the shortest hyperpath problem we consider here: Given a cell-signaling network whose reactants and reactions are modeled by the vertices and weighted hyperedges of a directed hypergraph, together with a set of sources and a target, find a hyperpath consisting of hyperedges from the sources to the target of minimum total weight.

## Installation

Installation requires the cplex and halp packages. Building the hypergraph files requires a biopax parser from Ritz et al. available here: https://github.com/annaritz/pathway-connectivity.

## Running the code

To run the code, use the following on the command line:
`python run.py --name <dataset name> --cyclicparallel --target <target name>`
where dataset name is the name given to the -hyperedges.txt and -hypernodes.txt files,
and target name is the name of the target in the same files.
