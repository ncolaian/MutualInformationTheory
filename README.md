# Mutual Information Theory Analysis
Code to perform a full MIT analysis


### Quick information about MIT
MIT relies on Shannon Entropy, or the measurement of randomness. Two residues in a multiple sequence alignmnet (MSA) can be identified as "co-evolving" by measuring how much randomness in a particulare residue of a protein can be explained by another. Put simply, if the identity of a residue at position 3 can be used to predict the residue at position 5, 7, etc. 

This type of measurement is useful to calculate within residue contacts and interaction sites. I have used it in the past on flagellin genes while studying two apposing evoltionary forces - motility vs detection by the plant immune system. By studying bacteria coming from plants, we found co-evolving residues that allow flagellin to retain motility while reducing detection by the plant immune system. 

This was really cool, and what initially got me interested in the algorithm. However, I started becoming more interested in using the algorith for understanding protein interactions. In the plant immune system the recognition of effectors or defense elliciting molecules requires the coordinated interaction of multiple proteins. I wondered if this algorithm could not only identify the regions of the protein that interact with the defense elliciting molecule but also its response partner. Part of this story was figured out by _ et al. where they used entropy to identify the binding sites of the protein and then classify the different response proteins. However, one big mystery remains for many defense proteins - who are their interacting partners and what molecule allows them to physically interact?

In this end, I was hoping to use mutual information theory as a part of analysis to predict if defense proteins interact with each other. My rationale is that changes in one protein to respond to evolution of its defense elliciting molecule would likely need changes in its binding partner, thus creating co-evolution. I figured that proteins that have more significant co-evolving residues than whould be seen by random chance would likely indicate interactions between proteins. This is significant because it is extremely hard to test protein interaction networks in-vivo like shown by Mott et al. (HERE) and you generally need to identify the elliciting molecule to identify a true relationship, however if you have information about potential partnerships, you may be able to bait out the required binding molecule.

## Requirements before analysis
In order to perform the analysis, one needs to identify homologs of your protein of interest across 100's-1000's of genomes. 

If you want to perform MIT on a single protein, this is good enough. If you want to test coevolution across proteins you need to make sure that the MSA's are ordered by genome. If there are multiple copies of a single gene you can either repeat the single copy gene or remove the genome from the analysis. I would lean towards removal only because with multiple copies it can allow one homolog to be non-functional. If you want to be more thorough, you can use gene synteny to determine the true functional gene. 

When graphing the data, one normally uses a heatmap to pull the data together. If you want to make the graph meaningful, you should renumber the sequence based on a well known reference sequence. This makes interpretation of the graph much easier - especially if you can use Alphafold to predict the conformational shape of the protein or even better conformational data based on crystallography. 
