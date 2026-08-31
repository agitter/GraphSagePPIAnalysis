# Paper method and URL extracts

## OhmNet.txt

### URLs
- `http://creativecommons.org/licenses/by-nc/4.0/),`
- `http://snap.stanford.`
### Gene Ontology
```
 Leskovec, 2016): This approach learns d-
work of the global PPI network, induced by the set of specifically
co-expressed edges in that tissue.                                               dimensional features for nodes based on a biased random walk
                                                                                 procedure that flexibly explores network neighborhoods of nodes.

4.3 Tissue-specific cellular functions and gene                              In addition, we evaluate the performance of OhmNet against the
                                                                             following tissue-specific/agnostic function prediction approaches:
annotations
Associations between tissues and cellular functions were retrieved           •   GeneMania (Zuberi et al., 2013): This is a supervised approach
from Greene et al. (2015). Greene et al. manually curated biological             that takes a multi-layer network as input and directly predicts
processes in the Gene Ontology (GO) (Ashburner et al., 2000) and                 cellular functions in two separate phases. In the first phase, it ag-
mapped them to tissues in the BRENDA Tissue Ontology (Chang                      gregates the layers into one weighted network by weighting the
et al., 2014) based on whether a given biological process is specific-           layers according to their utility for predicting a given function. It
ally active in a given tissue. The data is provided as a supplementar            then uses a label propagation algorithm on the weighted network
dataset in Greene et al. (2015). An example of a cellular function-              to predict the function.
tissue pair is "low-density lipoprotein particle remodeling" in the          •   Tissue-specific network propagation (Magger et al., 2012): This ap-
blood plasma tissue.                                                             proach assigns a prior score to proteins associated with known func-
    All gene annotations were propagated along the ontology hier-                tions that are phenotypically similar to the query function. This score
archy. Considered are functions with at least 15 annotated proteins              is then propagated through a network in an iterative process. The ap-
(Guan et al., 2012). In total, there are 584 tissue-specific cellular            proach was developed for tissue-specific disease gene prioritization.
functions covering 48 distinct tissues. Each tissue-specific func-           •   Network-based tissue-specific support vector machine (SVM)
tion is assigned to one or more leaves in the tissue hierarchy                   (Guan et al., 2012): This approach adopts the network-based
(Section 4.1).                                                                   candidate gene prediction scheme. Essentially, the connection
                                                                                 weights in a network to all positive examples (i.e. genes already
                                                                                 known to be related to a phenotype) are utilized as features for
5 Results                                                                        linear SVM classification. The approach was developed for
The OhmNet’s objective in Equation (5) is independent of any                     tissue-specific phenotype and disease gene prioritization.
downstream task. This flexibility offered by OhmNet makes the
                     
```
```
archy. Again, OhmNet                    Conflict of Interest: none declared.
i198                                                                                                                                    M.Zitnik and J.Leskovec


References                                                                             Magger,O. et al. (2012) Enhancing the prioritization of disease-causing genes
                                                                                          through tissue specific protein interaction networks. PLoS Comput. Biol., 8,
Antanaviciute,A. et al. (2015) GeneTIER: prioritization of candidate disease              e1002690.
  genes using tissue-specific gene expression profiles. Bioinformatics, 31,            Menche,J. et al. (2015) Uncovering disease-disease relationships through the
  2728–2735.                                                                              incomplete interactome. Science, 347, 1257601.
Ashburner,M. et al. (2000) Gene Ontology: tool for the unification of biology.         Mikolov,T. et al. (2013) Efficient estimation of word representations in vector
  Nat. Genet., 25, 25–29.                                                                 space. arXiv:1301.3781.
Barutcuoglu,Z. et al. (2006) Hierarchical multi-label prediction of gene func-         Mostafavi,S., and Morris,Q. (2009) Using the gene ontology hierarchy when
  tion. Bioinformatics, 22, 830–836.                                                      predicting gene function. In UAI, AUAI Press, Corvallis, pp. 419–427.
Belkin,M., and Niyogi,P. (2001) Laplacian eigenmaps and spectral techniques            Mostafavi,S. et al. (2008) GeneMANIA: a real-time multiple association network
  for embedding and clustering. In NIPS, vol. 14, MIT Press, Cambridge,                   integration algorithm for predicting gene function. Genome Biol., 9, 1.
  pp. 585–591.                                                                         Nickel,M. et al. (2011) A three-way model for collective learning on multi-
Cannistraci,C.V. et al. (2013) Minimum curvilinearity to enhance topological              relational data. In ICML, ACM, Bellevue, pp. 809–816.
  prediction of protein interactions by network embedding. Bioinformatics,             Okabe,Y., and Medzhitov,R. (2014) Tissue-specific signals control reversible
  29, i199–i209.                                                                          program of localization and functional polarization of macrophages. Cell,
Carvunis,A.-R., and Ideker,T. (2014) Siri of the cell: what biology could learn           157, 832–844.
  from the iPhone. Cell, 157, 534–538.                                                 Orchard,S. et al. (2013) The MIntAct projectintact as a common curation plat-




                                                                                                                                                                          Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
Chang,A. et al. (2014) BRENDA in 2015: exciting developments in its 25th                  form for 11 molecular interaction databases., Nucleic Acids Res. 42,
  year of existence. Nucleic Acids Res. 43, D439–D446.                                    D358–D363.
Chatr-Aryamontri,A. et al. (2015) The BioGRID interaction database: 2015               Perozzi,B. et al. (2014) Deepwalk: online learning of social representations. In
  
```
```
ease-causing genes
                                                                                          through tissue specific protein interaction networks. PLoS Comput. Biol., 8,
Antanaviciute,A. et al. (2015) GeneTIER: prioritization of candidate disease              e1002690.
  genes using tissue-specific gene expression profiles. Bioinformatics, 31,            Menche,J. et al. (2015) Uncovering disease-disease relationships through the
  2728–2735.                                                                              incomplete interactome. Science, 347, 1257601.
Ashburner,M. et al. (2000) Gene Ontology: tool for the unification of biology.         Mikolov,T. et al. (2013) Efficient estimation of word representations in vector
  Nat. Genet., 25, 25–29.                                                                 space. arXiv:1301.3781.
Barutcuoglu,Z. et al. (2006) Hierarchical multi-label prediction of gene func-         Mostafavi,S., and Morris,Q. (2009) Using the gene ontology hierarchy when
  tion. Bioinformatics, 22, 830–836.                                                      predicting gene function. In UAI, AUAI Press, Corvallis, pp. 419–427.
Belkin,M., and Niyogi,P. (2001) Laplacian eigenmaps and spectral techniques            Mostafavi,S. et al. (2008) GeneMANIA: a real-time multiple association network
  for embedding and clustering. In NIPS, vol. 14, MIT Press, Cambridge,                   integration algorithm for predicting gene function. Genome Biol., 9, 1.
  pp. 585–591.                                                                         Nickel,M. et al. (2011) A three-way model for collective learning on multi-
Cannistraci,C.V. et al. (2013) Minimum curvilinearity to enhance topological              relational data. In ICML, ACM, Bellevue, pp. 809–816.
  prediction of protein interactions by network embedding. Bioinformatics,             Okabe,Y., and Medzhitov,R. (2014) Tissue-specific signals control reversible
  29, i199–i209.                                                                          program of localization and functional polarization of macrophages. Cell,
Carvunis,A.-R., and Ideker,T. (2014) Siri of the cell: what biology could learn           157, 832–844.
  from the iPhone. Cell, 157, 534–538.                                                 Orchard,S. et al. (2013) The MIntAct projectintact as a common curation plat-




                                                                                                                                                                          Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
Chang,A. et al. (2014) BRENDA in 2015: exciting developments in its 25th                  form for 11 molecular interaction databases., Nucleic Acids Res. 42,
  year of existence. Nucleic Acids Res. 43, D439–D446.                                    D358–D363.
Chatr-Aryamontri,A. et al. (2015) The BioGRID interaction database: 2015               Perozzi,B. et al. (2014) Deepwalk: online learning of social representations. In
  update. Nucleic Acids Res., 43, D470–D478.                                              KDD, ACM, pp. 701–710.
Costanzo,M. et al. (2016) A global genetic interaction network maps a wiring           Prasad,T.K. et al. (2009) Human protein reference database-2009 update.
  diagram of cellular function. Science, 353, aaf1420.                                    Nucleic Acids Res. 37, D76
```
```
          Nucleic Acids Res. 37, D767–D772.
De Domenico,M. et al. (2014) Navigability of interconnected networks under             Przulj,N. (2007) Biological network comparison using graphlet degree distri-
  random failures. PNAS, 111, 8351–8356.                                                  bution. Bioinformatics, 23, e177–e183.
De Domenico,M. et al. (2015) Ranking in interconnected multilayer networks             Radivojac,P. et al. (2013) A large-scale evaluation of computational protein
  reveals versatile nodes. Nat. Commun., 6, 6868.                                         function prediction. Nat. Methods, 10, 221–227.
De Domenico,M. et al. (2016) The physics of spreading processes in multilayer          Rakyan,V.K. et al. (2008) An integrated resource for genome-wide identifica-
  networks. Nat. Phys., 12, 901–906.                                                      tion and analysis of human tissue-specific differentially methylated regions
Dutkowski,J. et al. (2012) A gene ontology inferred from molecular networks.              (tdmrs). Genome Res., 18, 1518–1529.
  Nat. Biotechnol., 31, 38–45.                                                         Rolland,T. et al. (2014) A proteome-scale map of the human interactome net-
Fagerberg,L. et al. (2014) Analysis of the human tissue-specific expression by            work. Cell, 159, 1212–1226.
  genome-wide integration of transcriptomics and antibody-based prote-                 Ruepp,A. et al. (2010) CORUM: the comprehensive resource of mammalian
  omics. Mol. Cell. Proteom., 13, 397–406.                                                protein complexes-2009. Nucleic Acids Res. 38, D497–D501.
Ganegoda,G.U. et al. (2014) Prediction of disease genes using tissue-specified         Stojanova,D. et al. (2013) Using PPI network autocorrelation in hierarchical
  gene-gene network. BMC Syst. Biol., 8, S3.                                              multi-label classification trees for gene function prediction. BMC
Greene,C.S. et al. (2015) Understanding multicellular function and disease                Bioinformatics, 14, 1.
  with human tissue-specific networks. Nat. Genet., 47, 569–576.                       Tang,J. et al. (2015) Line: Large-scale information network embedding. In
Grover,A., and Leskovec,J. (2016) Node2vec: scalable feature learning for net-            WWW, pp. 1067–1077.
  works. In KDD, pp. 855–864.                                                          Tang,L. et al. (2012) Scalable learning of collective behavior. IEEE Trans.
GTEx,C. et al. (2015) The genotype-tissue expression (GTEx) pilot analysis:               Knowl. Data Eng., 24, 1080–1091.
  multitissue gene regulation in humans. Science, 348, 648–660.                        Tenenbaum,J.B. et al. (2000) A global geometric framework for nonlinear
Guan,Y. et al. (2012) Tissue-specific functional networks for prioritizing                dimensionality reduction. Science, 290, 2319–2323.
  phenotype and disease genes. PLoS Comput. Biol., 8, e1002694.                        Vidulin,V. et al. (2016) Extensive complementarity between gene function pre-
Hayes,W. et al. (2013) Graphlet-based measures are suitable for biological                diction methods. Bioinformatics, 32, 3645–3653.
  network comparison. Bioinformatics, 29, 483–491.                                     Wang,D. et al. (2016a) Structural deep network embedding. In KDD, ACM,
Hou,C. et al. (2014) Joint embedding learning and sparse regression: a       
```
### download
```
                                                                                                         Bioinformatics, 33, 2017, i190–i198
                                                                                                         doi: 10.1093/bioinformatics/btx252
                                                                                                                           ISMB/ECCB 2017




Predicting multicellular function through
multi-layer tissue networks
Marinka Zitnik and Jure Leskovec*
Department of Computer Science, Stanford University, Stanford, CA 94305, USA
*To whom correspondence should be addressed.




                                                                                                                                                                          Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
Abstract
Motivation: Understanding functions of proteins in specific human tissues is essential for insights
into disease diagnostics and therapeutics, yet prediction of tissue-specific cellular function remains
a critical challenge for biomedicine.
Results: Here, we present OhmNet, a hierarchy-aware unsupervised node feature learning ap-
proach for multi-layer networks. We build a multi-layer network, where each layer represents mo-
lecular interactions in a different human tissue. OhmNet then automatically learns a mapping of
proteins, represented as nodes, to a neural embedding-based low-dimensional space of features.
OhmNet encourages sharing of similar features among proteins with similar network neighbor-
hoods and among proteins activated in similar tissues. The algorithm generalizes prior work, which
generally ignores relationships between tissues, by modeling tissue organization with a rich multi-
scale tissue hierarchy. We use OhmNet to study multicellular function in a multi-layer protein inter-
action network of 107 human tissues. In 48 tissues with known tissue-specific cellular functions,
OhmNet provides more accurate predictions of cellular function than alternative approaches, and
also generates more accurate hypotheses about tissue-specific protein actions. We show that tak-
ing into account the tissue hierarchy leads to improved predictive power. Remarkably, we also
demonstrate that it is possible to leverage the tissue hierarchy in order to effectively transfer cellu-
lar functions to a functionally uncharacterized tissue. Overall, OhmNet moves from flat networks to
multiscale models able to predict a range of phenotypes spanning cellular subsystems.
Availability and implementation: Source code and datasets are available at http://snap.stanford.
edu/ohmnet.
Contact: jure@cs.stanford.edu



1 Introduction                                                                        functional information from protein interaction networks lack tissue
A unified view of human diseases and cellular functions across a                      specificity as they assume that cellular function is constant across
broad range of human tissues is essential, not only for understanding                 organs and tissues (Barutcuoglu et al., 2006; Kramer et al., 2014;
basic biology but also for interpreting genetic variation and develop-                Mostafavi et al.,
```
```
disease associations. As such, these approaches account for tissue                 resentations for proteins that are consistent with the tissue
specificity, but they do not resolve the challenge of predicting gene–             hierarchy.
function relationships that might be specific to a particular tissue.                  Our experiments focus on three tasks defined on a multi-layer
To be able to predict a range of tissue-specific functions one needs               tissue network: (1) a multi-label node classification task, where
to design scalable multiscale models that can relate tissues to each               every protein is assigned zero, one or more tissue-specific cellular
other, extract rich feature representations for proteins in each tissue-           functions; (ii) a transfer learning task, where we predict cellular




                                                                                                                                                            Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
specific network, and then use the extracted features for tissue-                  functions for a protein in one tissue based on classifiers trained on
specific cellular function prediction.                                             features from other tissues; and (iii) a network-embedding visualiza-
                                                                                   tion task, where we create meaningful tissue-specific visualizations
                                                                                   that lay out proteins on a 2D space. Since the multiscale protein fea-
1.1 Present work                                                                   ture vectors returned by OhmNet are task-independent, we use
We present OhmNet, an algorithm for hierarchy-aware unsuper-                       OhmNet one time only to learn the features for proteins in every tis-
vised feature learning in multi-layer networks. Our focus is on learn-             sue and at every scale of the tissue hierarchy. We can then solve the
ing features of proteins in different tissues. We represent each tissue            cellular function prediction task for any tissue using the appropriate
as a network, where nodes represent proteins. Tissue networks act                  tissue-specific protein features.
as layers in a multi-layer network, where we use a hierarchy to                        We contrast OhmNet’s performance with that of state-of-the-art
model dependencies between the layers (i.e. tissues) (Fig. 1). We                  approaches for feature learning (Cannistraci et al., 2013; Grover
then develop a computational framework that learns features of                     and Leskovec, 2016; Nickel et al., 2011; Tang et al., 2015;),
each node (i.e. protein) by taking into consideration connections be-              approaches for tissue-independent cellular function prediction
tween the nodes within each layer, together with inter-layer relation-             (Mostafavi et al., 2008; Zuberi et al., 2013), and approaches for pri-
ships between proteins active on different layers. More precisely,                 oritization of disease-causing genes in tissue-specific protein inter-
our approach embeds each protein in each tissue in a d-dimensional                 action networks (Guan et al., 2012; Magger et al., 2012), which we
feature space such that p
```
```
y directly optimize the objective function for a downstream
et al., 2012), or nonlinear techniques based on multi-dimensional         prediction task, such as cellular function prediction in a particular
scaling (Belkin and Niyogi, 2001; Hou et al., 2014; Tenenbaum             tissue, using several layers of nonlinear transformations. Second,
et al., 2000). These methods have two important drawbacks. First,         those architectures do not model rich graph structures, such as
they do not account for important structures typically exhibited in       multi-layer networks with hierarchies.
networks such as high sparsity and skewed degree distribution.
Second, matrix factorization methods perform a global factorization
of the data matrix while a local-centric method might often yield         3 Feature learning in multi-layer networks




                                                                                                                                                        Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
more useful feature representations (Kramer et al., 2014).                We formulate feature learning in multi-layer networks as a max-
     Limitations of matrix factorization are overcome by neural net-      imum likelihood optimization problem. Let V be a given set of N
work embeddings. Recent studies focused on embedding nodes into           nodes (e.g. proteins) fu1 ; u2 ; . . . ; uN g; and let there be K types of
low-dimensional vector spaces by first using random walks to con-         edges (e.g. protein interactions in different tissues) between pairs of
struct the network neighborhood of every node in the graph, and           nodes u1 ; u2 ; . . . ; uN . A multi-layer network is a general system in
then optimizing an objective function with network neighborhoods          which each biological context is represented by a distinct layer i
as input (Grover and Leskovec, 2016; Perozzi et al., 2014; Tang           (where i ¼ 1; 2; . . . ; K) of a system (Fig. 1). We use the term single-
et al., 2015). The objective function is carefully designed to preserve   layer network (layer) for the network Gi ¼ ðVi ; Ei Þ that indicates the
both the local and global network structures. A state-of-the-art neu-     edges Ei between nodes Vi  V within the same layer i. Our analysis
ral network embedding algorithm is the Node2vec algorithm                 is general and applies to any (un)directed, (un)weighted multi-layer
(Grover and Leskovec, 2016), which learns feature representations         network.
as follows: it scans over the nodes in a network, and for every node          We take into account the possibility that a node uk from layer i
it aims to embed it such that the node’s features can predict nearby      can be related to any other node uh in any other layer j. We encode
nodes, that is, node’s feature predict which other nodes are part of      information about the dependencies between layers in a hierarchical
its network neighborhood. Node2vec can explore different network          manner that we use in the learning process. Let the hierarchy be a
neighborhoods to embed nodes based on the principles of homo-             directed tree M defined over a set M of elements by the parent-child
phily (i.e. network communities) as well as structural equivalence        relationships given by p : M ! M; where pðiÞ is the parent of elem-
(i.e. structural roles of node
```
```
archy elements toward fea-
node pair is modeled independently as PrðNi ðuÞjfi ðuÞÞ ¼                    tures in the common parent element in the hierarchy.
Q
   v2Ni ðuÞ Prðvjfi ðuÞÞ. The conditional likelihood is a softmax unit
parameterized by a dot product of nodes’ features, which corres-
                                                                             3.2.1 Node features at multiple scales
ponds to a single-layer feed-forward neural network: Prðvjfi ðuÞÞ ¼
                     P                                                       It is important to notice that OhmNet’s structured regularization
exp ðfi ðvÞfi ðuÞÞ= z2Vi exp ðfi ðzÞfi ðuÞÞ. Given a node u, maximiza-
                                                                             allows us to learn feature representations at multiple scales. For ex-




                                                                                                                                                           Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
tion of xi ðuÞ tries to maximize classification of nodes in u’s network
                                                                             ample, consider a multi-layer network in Figure 2, consisting of four
neighborhood based on u’s learned representation.
                                                                             layers that are interrelated by a two-level hierarchy. OhmNet learns
   The objective Xi is defined for each layer i:                             the mappings fi, fj, fk and fl that map nodes in each layer into a d-di-
                     X                                                       mensional feature space. In addition, OhmNet also learns the map-
              Xi ¼       xi ðuÞ; for i ¼ 1; 2; . . . ; K:             (2)
                      u2Vi
                                                                             ping f2 representing features for nodes appearing in the hierarchy
                                                                             leaves T2, i.e. Vi [ Vj , at an intermediate scale, and the mapping f1
The objective is inspired by the intuition that nodes with similar net-      representing features for nodes appearing in the hierarchy leaves T1,
work neighborhoods tend to have similar meanings, or roles, in a             i.e. Vi [ Vj [ Vk [ Vl , at the highest scale.
network. It formalizes this intuition by encouraging nodes in similar             The modeling of relationships between layers in a multi-layer
network neighborhoods to share similar features.                             network has several implications:
    We found that a flexible notion of a network neighborhood Ni is
                                                                             •   First, the model encourages nodes which are in nearby layers in
crucial to achieve excellent predictive accuracy on a downstream
cellular function prediction task (Grover and Leskovec, 2016). For               the hierarchy to share similar features.
                                                                             •   Second, the model shares statistical strength across the hierarchy
that reason, we use a randomized procedure to sample many differ-
ent neighborhoods of a given node u. Technically, the network                    as nodes in different layers representing the same protein share
neighborhood Ni ðuÞ
```
```
          f1 ;f2 ;...;fjMj
                                           i2T   j2M
                                                                              OhmNet algorithm scales to large multi-layer networks because
which includes the single-layer network objectives for all network            each phase is parallelizable and executed asynchronously. The
layers, and the hierarchical dependency objectives for all hierarchy          choice to model the dependencies between network layers using the
elements. In Equation (5), parameter k is a user-specified parameter          hierarchical model requires OðjMjNÞ time instead of the fully pair-
representing the regularization strength. While the optimization              wise model, which requires OðK2 NÞ time.
problem in Equation (5) is non-convex due to the non-convexity of




                                                                                                                                                                 Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
the single-layer objective (Grover and Leskovec, 2016), stochastic
gradient with negative sampling can be used to efficiently solve the          4 Tissue-specific interactome data
problem.
                                                                              To construct the human PPI network, tissue-specific network layers,
       One appealing property of OhmNet is that by solving the prob-
                                                                              tissue hierarchy and tissue-specific gene–function relationships, we
lem in Equation (5) we obtain estimates for functions f1 ; f2 ; . . . ; fK
                                                                              downloaded and used standard protein, tissue and function informa-
located in the leaf elements of the hierarchy (i.e. layers of a given
                                                                              tion from various reputable data sources.
multi-layer network), as well as estimates for functions fKþ1 ; fKþ2 ;
. . . ; fjMj located in the internal elements of the hierarchy.
                                                                              4.1 Tissue hierarchy
                                                                              We retrieved the mapping of tissues in the Human Protein Reference
3.4 The OhmNet algorithm
                                                                              Database (HPRD) (Prasad et al., 2009) to tissues in the BRENDA
The pseudocode for OhmNet is given in Algorithm 1.
                                                                              Tissue Ontology (Chang et al., 2014) from Greene et al. (2015). The
    In the first phase, OhmNet applies the Node2vec’s algorithm
                                                                              data is provided as a supplementary dataset in Greene et al. (2015).
(Grover and Leskovec, 2016) to construct network neighborhoods
                                                                              The hierarchical relationships between tissues were then determined
for each node in every layer. Given a layer Gi and a node u 2 Vi , the
                                                                              by the directed acyclic graph structure of the BRENDA Tissue
algorithm simulates a user-defined number of fixed length random
        
```
```
 due to the non-convexity of




                                                                                                                                                                 Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
the single-layer objective (Grover and Leskovec, 2016), stochastic
gradient with negative sampling can be used to efficiently solve the          4 Tissue-specific interactome data
problem.
                                                                              To construct the human PPI network, tissue-specific network layers,
       One appealing property of OhmNet is that by solving the prob-
                                                                              tissue hierarchy and tissue-specific gene–function relationships, we
lem in Equation (5) we obtain estimates for functions f1 ; f2 ; . . . ; fK
                                                                              downloaded and used standard protein, tissue and function informa-
located in the leaf elements of the hierarchy (i.e. layers of a given
                                                                              tion from various reputable data sources.
multi-layer network), as well as estimates for functions fKþ1 ; fKþ2 ;
. . . ; fjMj located in the internal elements of the hierarchy.
                                                                              4.1 Tissue hierarchy
                                                                              We retrieved the mapping of tissues in the Human Protein Reference
3.4 The OhmNet algorithm
                                                                              Database (HPRD) (Prasad et al., 2009) to tissues in the BRENDA
The pseudocode for OhmNet is given in Algorithm 1.
                                                                              Tissue Ontology (Chang et al., 2014) from Greene et al. (2015). The
    In the first phase, OhmNet applies the Node2vec’s algorithm
                                                                              data is provided as a supplementary dataset in Greene et al. (2015).
(Grover and Leskovec, 2016) to construct network neighborhoods
                                                                              The hierarchical relationships between tissues were then determined
for each node in every layer. Given a layer Gi and a node u 2 Vi , the
                                                                              by the directed acyclic graph structure of the BRENDA Tissue
algorithm simulates a user-defined number of fixed length random
                                                                              Ontology. Examples of tissues included: muscle, adrenal cortex,
walks started at node u (step 4 in Algorithm 1).
                                                                              bone marrow and spleen (Fig. 3).
    In the second phase, OhmNet uses an iterative approach in
which features associated with each object in the hierarchy are itera-
tively updated by fixing the rest of the features. The iterative              4.2 Tissue-specific interaction networks
                                                                              We took the gene-to-tissue mapping compiled by Greene et al.
                                                                              (2015). Greene et al. mapped genes to HPRD 
```
```
                                                                        for protein interaction prediction, aiming to embed protein pairs
every edge in the global PPI network was labeled as specifically co-
                                                                                 representing good candidate interactions closer to each other. It
expressed in that tissue using the criterion developed by Greene
et al. (2015). Greene et al. labeled each edge as specifically co-               utilizes a network denoising method as well as structural infor-
expressed if either both proteins are specific to that tissue or one             mation provided by the PPI network topology.
                                                                             •   LINE (Tang et al., 2015): This approach first learns d=2 dimensions




                                                                                                                                                            Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
protein is tissue-specific and the other is ubiquitous. Lists of specific-
ally co-expressed proteins were retrieved from Greene et al. (2015).             based on immediate network neighbors of nodes, and then the next
Finally, the PPI network specific to a particular tissue is a subnet-            d=2 dimensions based on network neighbors at a 2-hop distance.
                                                                             •   Node2vec (Grover and Leskovec, 2016): This approach learns d-
work of the global PPI network, induced by the set of specifically
co-expressed edges in that tissue.                                               dimensional features for nodes based on a biased random walk
                                                                                 procedure that flexibly explores network neighborhoods of nodes.

4.3 Tissue-specific cellular functions and gene                              In addition, we evaluate the performance of OhmNet against the
                                                                             following tissue-specific/agnostic function prediction approaches:
annotations
Associations between tissues and cellular functions were retrieved           •   GeneMania (Zuberi et al., 2013): This is a supervised approach
from Greene et al. (2015). Greene et al. manually curated biological             that takes a multi-layer network as input and directly predicts
processes in the Gene Ontology (GO) (Ashburner et al., 2000) and                 cellular functions in two separate phases. In the first phase, it ag-
mapped them to tissues in the BRENDA Tissue Ontology (Chang                      gregates the layers into one weighted network by weighting the
et al., 2014) based on whether a given biological process is specific-           layers according to their utility for predicting a given function. It
ally active in a given tissue. The data is provided as a supplementar            then uses a label propagation algorithm on the weighted network
dataset in Greene et al. (2015). An example of a cellular function-              to predict the function.
tissue pair is "low-density lipoprotein particle remodeling" in the          •   Tissue-specific network propagation (Magger et al., 2012): This ap-
blood plasma tissue.                                                             proach assigns a prior score t
```
```
         0.785 (60.030)             0.749 (60.032)
GeneMania                                 0.683 (60.077)     0.274 (60.094)      Hematopoietic stem cell           0.784 (60.035)             0.744 (60.036)
Network-based tissue-specific SVM         0.701 (60.091)     0.281 (60.059)      Blood plasma                      0.784 (60.027)             0.703 (60.039)
Tissue-specific network propagation       0.675 (60.051)     0.265 (60.083)      Smooth muscle                     0.778 (60.031)             0.729 (60.041)
OhmNet (Section 3)                        0.756 (60.067)     0.336 (60.045)      Average                           0.799                      0.746

   Values in the brackets are halves of the interquartile distance. OhmNet’s        Shown are the scores for ten tissues with best performance on cellular func-




                                                                                                                                                                    Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
results are statistically significant with a P-value of < 0.05                   tion prediction task. “Non-transfer”: a classifier is trained on a target tissue
                                                                                 and then used to predict cellular functions in the same tissue (Section 5.1).
                                                                                 “Transfer”: classifiers are trained on all non-target tissues and then used to
                                                                                 predict cellular functions in the target tissue (Section 5.2).


                                                                                 generate more accurate hypotheses about tissue-specific protein
                                                                                 actions.


                                                                                 5.2 Transfer of cellular functions to a new tissue
                                                                                 5.2.1 Experimental setup
                                                                                 In the transfer learning setting, we attempt to transfer knowledge
                                                                                 learned in one or more source layers and use it for prediction in a
                                                                                 target layer.
                                                                                     As before, we apply OhmNet to obtain a separate feature vector
                                                                                 for every node and every layer in an unsupervised way. We then con-
                                                                                 sider, in turn, every tissue as a target layer and all other tissues as
                                                                                 source layers. For every function and every source layer, we train a
                                                                                 separate classifier using the same classification model as in Section
Fig. 4. Area under ROC curve (AUROC) scores for tissue-specific cellular func-   5.1. We then predict functions for the target layer using only classi-
tion pred
```
```
                to 20.3% over the closest benchmark in AUC scores (scores not
approaches able to directly profile proteins’ distinct interaction               shown). Notice that we exclude GeneMania in the comparison be-
neighborhoods in different tissues can leverage this specificity to              cause it is not amenable to transfer learning. This result suggests
Predicting multicellular function through multi-layer tissue networks                                                                                           i197


   A                                                        B                                                        C




Fig. 5. Visualization of the brain tissue-specific protein interaction networks. (A) The two-level brain tissue hierarchy as specified by the BRENDA Tissue Ontology




                                                                                                                                                                         Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
(Chang et al., 2014) and used in the case study in Section 5.3. Leaves of the hierarchy (in blue) represent nine brain tissues each of which is associated with a tis-
sue-specific protein interaction network. (B) Visualization of the brainstem-specific networks. The proteins are mapped to the 2D space using the t-SNE package
with learned features as input. Color of a node indicates the tissue of the protein. (C) Visualization of the brain-specific networks. The proteins are mapped and
colored using the same procedure as in B


that considering the relationships between tissues when learning                      produces a meaningful layout of the nine brain tissue-specific
features for proteins has a significant impact on transfer performance.               networks.
    In general, we observed that the transferability of classifiers                       In addition, we repeated this analysis by visualizing protein fea-
decreased when the tree-based distance between the source and the                     tures learned by running principal component analysis (PCA) or
target tissue in the tissue hierarchy increased, which is consistent with             non-negative matrix factorization (NMF) algorithm on the brain-
the empirical evidence in transfer learning (Yosinski et al., 2014). This             specific PPI networks. Acknowledging the subjective nature of this
also matches our intuition that a source tissue should be most inform-                analysis, we observed that visualizations using PCA or NMF were
ative for predicting cellular functions in an anatomically close target               not very meaningful, as proteins belonging to the same tissue were
tissue (e.g. source and target tissues are both part of the same organ).              not clustered together (data not shown).
                                                                                          OhmNet’s result in Figure 5 is especially appealing because of
5.3 The multiscale model of brain tissues                                             two reasons. First, it shows that OhmNet can learn node features
We have seen in Section 4.1 that human tissues have a multi-level                     that adhere to a given hierarchy of layers. In the brain example,
hierarchical organization. The tissue hierarchy categorizes tissues                   OhmNet learns the protein features that expose 
```
```
         Nickel,M. et al. (2011) A three-way model for collective learning on multi-
Cannistraci,C.V. et al. (2013) Minimum curvilinearity to enhance topological              relational data. In ICML, ACM, Bellevue, pp. 809–816.
  prediction of protein interactions by network embedding. Bioinformatics,             Okabe,Y., and Medzhitov,R. (2014) Tissue-specific signals control reversible
  29, i199–i209.                                                                          program of localization and functional polarization of macrophages. Cell,
Carvunis,A.-R., and Ideker,T. (2014) Siri of the cell: what biology could learn           157, 832–844.
  from the iPhone. Cell, 157, 534–538.                                                 Orchard,S. et al. (2013) The MIntAct projectintact as a common curation plat-




                                                                                                                                                                          Downloaded from academic.oup.com/bioinformatics/article/33/14/i190/3953967 by guest on 16 August 2026
Chang,A. et al. (2014) BRENDA in 2015: exciting developments in its 25th                  form for 11 molecular interaction databases., Nucleic Acids Res. 42,
  year of existence. Nucleic Acids Res. 43, D439–D446.                                    D358–D363.
Chatr-Aryamontri,A. et al. (2015) The BioGRID interaction database: 2015               Perozzi,B. et al. (2014) Deepwalk: online learning of social representations. In
  update. Nucleic Acids Res., 43, D470–D478.                                              KDD, ACM, pp. 701–710.
Costanzo,M. et al. (2016) A global genetic interaction network maps a wiring           Prasad,T.K. et al. (2009) Human protein reference database-2009 update.
  diagram of cellular function. Science, 353, aaf1420.                                    Nucleic Acids Res. 37, D767–D772.
De Domenico,M. et al. (2014) Navigability of interconnected networks under             Przulj,N. (2007) Biological network comparison using graphlet degree distri-
  random failures. PNAS, 111, 8351–8356.                                                  bution. Bioinformatics, 23, e177–e183.
De Domenico,M. et al. (2015) Ranking in interconnected multilayer networks             Radivojac,P. et al. (2013) A large-scale evaluation of computational protein
  reveals versatile nodes. Nat. Commun., 6, 6868.                                         function prediction. Nat. Methods, 10, 221–227.
De Domenico,M. et al. (2016) The physics of spreading processes in multilayer          Rakyan,V.K. et al. (2008) An integrated resource for genome-wide identifica-
  networks. Nat. Phys., 12, 901–906.                                                      tion and analysis of human tissue-specific differentially methylated regions
Dutkowski,J. et al. (2012) A gene ontology inferred from molecular networks.              (tdmrs). Genome Res., 18, 1518–1529.
  Nat. Biotechnol., 31, 38–45.                                                         Rolland,T. et al. (2014) A proteome-scale map of the human interactome net-
Fagerberg,L. et al. (2014) Analysis of the human tissue-specific expression by            work. Cell, 159, 1212–1226.
  genome-wide integration of transcriptomics and antibody-based prote-                 Ruepp,A. et al. (2010) CORUM: the comprehensive resource of mammalian
  omics. Mol. Cell. Proteom., 13, 397–406.                                                
```
## Greene2015_sup.txt

## Greene2015.txt

### URLs
- `http://giant.princeton.edu/;`
- `http://libsleipnir.bitbucket.org/;`
- `http://tribe.`
- `http://www.nature.com/`
### Sleipnir
```
eans to                        biogenesis. PLoS Comput. Biol. 5, e1000322 (2009).
                                                   address these challenges.                                                        13. Park, C.Y. et al. Functional knowledge transfer for high-accuracy prediction of
                                                                                                                                        under-studied biological processes. PLoS Comput. Biol. 9, e1002957 (2013).
                                                                                                                                    14. Jansen, R. et al. A Bayesian networks approach for predicting protein-protein
                                                   URLs. GIANT, a web portal for tissue-specific functional networks,                   interactions from genomic data. Science 302, 449–453 (2003).
                                                   http://giant.princeton.edu/; Sleipnir, an open source library for                15. Lee, I., Date, S.V., Adai, A.T. & Marcotte, E.M. A probabilistic functional network
                                                                                                                                        of yeast genes. Science 306, 1555–1558 (2004).
                                                   functional genomics, http://libsleipnir.bitbucket.org/; Tribe, a web             16. Mostafavi, S., Ray, D., Warde-Farley, D., Grouios, C. & Morris, Q. GeneMANIA: a
                                                   service that provides cross-server analysis of gene sets, http://tribe.              real-time multiple association network integration algorithm for predicting gene
                                                   greenelab.com/.                                                                      function. Genome Biol. 9 (suppl. 1), S4 (2008).
                                                                                                                                    17. Hwang, S., Rhee, S.Y., Marcotte, E.M. & Lee, I. Systematic prediction of gene
                                                                                                                                        function in Arabidopsis thaliana using a probabilistic functional gene network.
                                                   Methods                                                                              Nat. Protoc. 6, 1429–1442 (2011).
                                                   Methods and any associated references are available in the online                18. Kofler, S., Nickel, T. & Weis, M. Role of cytokines in cardiovascular diseases:
                                                                                                                                        a focus on endothelial responses to inflammation. Clin. Sci. 108, 205–213
                                                   version of the paper.                                                                (2005).
                                                                                                                                    19. Liu, J.Z. et al. A versatile gene-based test for genome-wide association studies.
                                                                                                                                        Am. J. Hum. Genet. 87, 139–145 (2010).
     
```
```
                             under-studied biological processes. PLoS Comput. Biol. 9, e1002957 (2013).
                                                                                                                                    14. Jansen, R. et al. A Bayesian networks approach for predicting protein-protein
                                                   URLs. GIANT, a web portal for tissue-specific functional networks,                   interactions from genomic data. Science 302, 449–453 (2003).
                                                   http://giant.princeton.edu/; Sleipnir, an open source library for                15. Lee, I., Date, S.V., Adai, A.T. & Marcotte, E.M. A probabilistic functional network
                                                                                                                                        of yeast genes. Science 306, 1555–1558 (2004).
                                                   functional genomics, http://libsleipnir.bitbucket.org/; Tribe, a web             16. Mostafavi, S., Ray, D., Warde-Farley, D., Grouios, C. & Morris, Q. GeneMANIA: a
                                                   service that provides cross-server analysis of gene sets, http://tribe.              real-time multiple association network integration algorithm for predicting gene
                                                   greenelab.com/.                                                                      function. Genome Biol. 9 (suppl. 1), S4 (2008).
                                                                                                                                    17. Hwang, S., Rhee, S.Y., Marcotte, E.M. & Lee, I. Systematic prediction of gene
                                                                                                                                        function in Arabidopsis thaliana using a probabilistic functional gene network.
                                                   Methods                                                                              Nat. Protoc. 6, 1429–1442 (2011).
                                                   Methods and any associated references are available in the online                18. Kofler, S., Nickel, T. & Weis, M. Role of cytokines in cardiovascular diseases:
                                                                                                                                        a focus on endothelial responses to inflammation. Clin. Sci. 108, 205–213
                                                   version of the paper.                                                                (2005).
                                                                                                                                    19. Liu, J.Z. et al. A versatile gene-based test for genome-wide association studies.
                                                                                                                                        Am. J. Hum. Genet. 87, 139–145 (2010).
                                                   Accession codes. Gene expression measurements of HASMCs with                     20. Keshava Prasad, T.S. et al. Human Protein Reference Database—2009 update.
                                                   and without IL-1β stimulation are available in the Gene Expression                   Nucleic Acids Res. 37, D767–D772 (2009).
                       
```
```
    mated the probability of global functional interactions for the tissue-naive



                                                   doi:10.1038/ng.3259                                                                                                                                 Nature Genetics
                                                   network. We assigned a prior probability of a functional relationship of 0.01       the connectivity measure described above. We evaluated the mean fold change
                                                   for all models, allowing edge probabilities to be compared across tissues.          of the top 20 returned results. We evaluated randomly selected matched size
                                                      Code availability. Integrations were performed with C++ naive Bayesian           sets of genes from each data set as controls.
                                                   learning implementations from the open source Sleipnir library for functional
                                                   genomics82.                                                                         Evaluation of tissue-specific processes, gene-level rewiring and disease-
                                                                                                                                       disease association. Mapping GO biological processes to tissues. To evaluate
                                                   Evaluation of tissue-specific functional relationships. We evaluated                tissue-specific functional rewiring in our networks, we needed associations
                                                   tissue-naive and tissue-specific functional networks using fivefold cross-          between tissues and tissue-specific processes. We used text matching followed
                                                   validation. The 6,062 genes represented in the tissue-specific knowledge-           by manual curation to map biological process (BP) terms in GO to tissue terms
                                                   base were randomly partitioned into 5 sets. For each cross-validation run,          in the BRENDA Tissue ontology (Supplementary Table 9).
                                                   gene pairs where neither gene was present in the holdout interval were used             Network connectivity of tissue-specific processes. For each tissue, we
                                                   for training. Any gene pair where both genes were present in the holdout            constructed a tissue-minus-naive network by subtracting edge probabilities
                                                   was used for evaluation of the AUC. The estimated performance of each of            of the naive network from those of the tissue network. Negative weights were
                                                   the 144 functional networks was summarized as the median AUC of the five            set to zero. In this subtracted network, positive scores corresponded to edges
                                                   cross-validation runs (Supplementary Table 8).                                      with a tissue network interaction probability greater than the naive network
                                                      Mapping data sets to tissues. We mapped data sets to tissues to compare with     probability. We expected relevant t
```
```
                  80. Burkard, T.R. et al. Initial characterization of the human central proteome. BMC
                                                   network specific to kidney, a tissue associated with hypertension54, where
                                                                                                                                               Syst. Biol. 5, 17 (2011).
                                                   the features of the classifier were the edge weights of the labeled examples to         81. Uhlen, M. et al. Towards a knowledge-based Human Protein Atlas. Nat. Biotechnol.
                                                   all the genes in the network. Genes were re-ranked using their distance                     28, 1248–1250 (2010).
                                                   from the hyperplane, which represented a network-based prioritization of                82. Huttenhower, C., Schroeder, M., Chikina, M.D. & Troyanskaya, O.G. The Sleipnir
                                                                                                                                               library for computational functional genomics. Bioinformatics 24, 1559–1561
                                                   a GWAS, termed NetWAS.
                                                                                                                                               (2008).
                                                      We applied NetWAS to a GWAS from the Women’s Genome Health Study                     83. Schmid, P.R., Palmer, N.P., Kohane, I.S. & Berger, B. Making sense out of massive
                                                   to identify additional genes involved in hypertension 51. The study focused                 data by going beyond differential expression. Proc. Natl. Acad. Sci. USA 109,
                                                   on three hypertension-related endpoints: systolic blood pressure, diastolic                 5594–5599 (2012).
                                                                                                                                           84. Aronson, A.R. Effective mapping of biomedical text to the UMLS Metathesaurus:
                                                   blood pressure and hypertension diagnosis. To calculate per-gene P values for
                                                                                                                                               the MetaMap program. Proc. AMIA Symp. 2001, 17–21 (2001).
npg




                                                   each endpoint, we used the versatile gene-based association study (VEGAS)               85. Bolstad, B.M., Irizarry, R.A., Astrand, M. & Speed, T.P. A comparison of normalization
                                                   system19. To generate a combined list across phenotypes, we combined results                methods for high density oligonucleotide array data based on variance and bias.
                                                   from each hypertension-related endpoint using summed ranks. Performance                     Bioinformatics 19, 185–193 (2003).
                                                                                                                                           86. Dai, M. et al. Evolving gene/transcript definitions significantly alter the interpretation
                        
```
### FKT
```
ther                  of 144 tissues (Supplementary Table 8) that contained at least ten C1 edges
                                                   condition were excluded from the set of negative examples and treated as                 between tissue-specific genes (T–T). This method incorporates the hierarchi-
                                                   neither related nor unrelated.                                                           cal relationships of tissues, allowing supervised methods to leverage these
                                                      Functional knowledge transfer. To increase the coverage of functional                 relationships.
                                                   interactions, we transferred experimentally confirmed mouse GO annota-                       Data integration. We constructed functional networks from genome-scale
                                                   tions to human functional analogs identified by FKT13, a high-specificity                data by performing a tissue-specific Bayesian integration. We trained one naive
                                                   annotation transfer method, for the 520 GO terms with mouse annotations.                 Bayesian classifier for each tissue using the tissue-specific standards described
                                                   This resulted in a tissue-naive gold standard of 604,038 functionally related            above. We also trained a classifier limited to only functional interactions to
                                                   gene pairs (positive examples) and 12,425,713 potentially unrelated pairs                generate a tissue-naive network. In each case, we constructed a class node, i.e.,
                                                   (negative examples).                                                                     the presence or absence of a functional relationship between a pair of genes
                                                      Ontology pruning. Gene-to-tissue annotations were obtained from the                   that is conditioned on nodes for each data set. For large-scale genomics data
                                                   Human Protein Reference Database (HPRD)20. HPRD tissues were mapped to                   sets, the assumption of conditional independence required for a naive Bayes
                                                   the BRENDA Tissue Ontology21 (BTO) using direct matching where possible                  classifier is often not met, so we calculated and corrected for non-biological
                                                   and manual curation where direct matches were unavailable (Supplementary                 conditional dependency13.
                                                   Table 7). Tissues with fewer than ten directly annotated genes were pruned                   Each tissue model trained on the hierarchy-aware tissue-specific knowledge
                                                   as non-informative from a molecular standpoint (for example, BTO:0001493,                was used to make genome-wide predictions by estimating the probability of
                                                   trunk). Pruning resulted in an ontology containing functional, as opposed                tissue-specific functional interaction between all pairs of genes. We also esti-
                             
```
### functional knowledge transfer
```
fect               (2003).
                                                   not just in the target tissue but also in all tissues. By disentangling          11. Myers, C.L. & Troyanskaya, O.G. Context-sensitive data integration and prediction
                                                   the functions of genes in specific tissues, integrated tissue-specific               of biological networks. Bioinformatics 23, 2322–2330 (2007).
                                                                                                                                    12. Hibbs, M.A. et al. Directing experimental biology: a case study in mitochondrial
                                                   networks learned from large data compendia present a means to                        biogenesis. PLoS Comput. Biol. 5, e1000322 (2009).
                                                   address these challenges.                                                        13. Park, C.Y. et al. Functional knowledge transfer for high-accuracy prediction of
                                                                                                                                        under-studied biological processes. PLoS Comput. Biol. 9, e1002957 (2013).
                                                                                                                                    14. Jansen, R. et al. A Bayesian networks approach for predicting protein-protein
                                                   URLs. GIANT, a web portal for tissue-specific functional networks,                   interactions from genomic data. Science 302, 449–453 (2003).
                                                   http://giant.princeton.edu/; Sleipnir, an open source library for                15. Lee, I., Date, S.V., Adai, A.T. & Marcotte, E.M. A probabilistic functional network
                                                                                                                                        of yeast genes. Science 306, 1555–1558 (2004).
                                                   functional genomics, http://libsleipnir.bitbucket.org/; Tribe, a web             16. Mostafavi, S., Ray, D., Warde-Farley, D., Grouios, C. & Morris, Q. GeneMANIA: a
                                                   service that provides cross-server analysis of gene sets, http://tribe.              real-time multiple association network integration algorithm for predicting gene
                                                   greenelab.com/.                                                                      function. Genome Biol. 9 (suppl. 1), S4 (2008).
                                                                                                                                    17. Hwang, S., Rhee, S.Y., Marcotte, E.M. & Lee, I. Systematic prediction of gene
                                                                                                                                        function in Arabidopsis thaliana using a probabilistic functional gene network.
                                                   Methods                                                                              Nat. Protoc. 6, 1429–1442 (2011).
                                                   Methods and any associated references are available in the online                18. Kofler, S., Nickel, T. & Weis, M. Role of cytokines in cardiovascular diseases:
         
```
```
                                    divided by the standard deviation. The resulting z scores were discretized                      U–T′].
                                                   into bins ((−infinity, −1.5), [–1.5, −0.5), [–0.5, 0.5), [0.5, 1.5), [1.5, 2.5), [2.5,     C3: negative functional edges between genes specifically coexpressed
                                                   3.5), [3.5, infinity)).                                                                        in the tissue [T–T and T–U].
                                                                                                                                              C4: negative functional edges between one gene expressed in the
                                                   Hierarchically aware knowledgebase construction via ontological pruning                        tissue and another specifically expressed in an unrelated tissue
                                                   with functional knowledge transfer.                                                            [T–T′ and U–T′].
                                                   Functional knowledge extraction. We constructed a tissue-naive functional
                                                   relationship gold standard from a set of 564 expert-selected GO biological               Among the four tissue classes, C1 represented tissue-specific functional
                                                   process terms and experimentally derived gene annotations (GO evidence                   relationships. To identify tissue-specific relationships, we constructed a
                                                   codes: EXP, IDA, IPI, IMP, IGI and IEP). Curators identified processes                   specific gold standard for each tissue by labeling edges in C1 as positives and
npg




                                                   testable through specific molecular experiments (Supplementary Table 6).                 edges in the other classes as negatives. Because C3 is defined on the basis of
                                                   Pairs of genes that were co-annotated to expert-selected terms after propaga-            tissue-expressed genes and C2 and C4 are defined on the basis of non-expressed
                                                   tion were treated as positive (i.e., functionally related) examples. Gene pairs          genes, the number of edges in these classes varied across tissues according
                                                   not co-annotated to any of these terms were considered as negative examples,             to how specific (cell type, tissue, organ or system), well studied (or easily
                                                   except in the following cases: (i) if two genes were annotated to two differ-            studied) and well curated (literature bias) they are. To construct comparable
                                                   ent GO terms with a significant number of shared genes (hypergeometric                   networks across tissues, we used a negative set composed of equal propor-
                                                   P value < 0.05) and (ii) if two genes were co-annotated to a set of ‘negative’           tions of edges from C2, C3 and C4. We limited all integrations to the set
                                                   GO terms that defined minimal relatedness77. Gene pairs that met either 
```
```
         networks across tissues, we used a negative set composed of equal propor-
                                                   P value < 0.05) and (ii) if two genes were co-annotated to a set of ‘negative’           tions of edges from C2, C3 and C4. We limited all integrations to the set
                                                   GO terms that defined minimal relatedness77. Gene pairs that met either                  of 144 tissues (Supplementary Table 8) that contained at least ten C1 edges
                                                   condition were excluded from the set of negative examples and treated as                 between tissue-specific genes (T–T). This method incorporates the hierarchi-
                                                   neither related nor unrelated.                                                           cal relationships of tissues, allowing supervised methods to leverage these
                                                      Functional knowledge transfer. To increase the coverage of functional                 relationships.
                                                   interactions, we transferred experimentally confirmed mouse GO annota-                       Data integration. We constructed functional networks from genome-scale
                                                   tions to human functional analogs identified by FKT13, a high-specificity                data by performing a tissue-specific Bayesian integration. We trained one naive
                                                   annotation transfer method, for the 520 GO terms with mouse annotations.                 Bayesian classifier for each tissue using the tissue-specific standards described
                                                   This resulted in a tissue-naive gold standard of 604,038 functionally related            above. We also trained a classifier limited to only functional interactions to
                                                   gene pairs (positive examples) and 12,425,713 potentially unrelated pairs                generate a tissue-naive network. In each case, we constructed a class node, i.e.,
                                                   (negative examples).                                                                     the presence or absence of a functional relationship between a pair of genes
                                                      Ontology pruning. Gene-to-tissue annotations were obtained from the                   that is conditioned on nodes for each data set. For large-scale genomics data
                                                   Human Protein Reference Database (HPRD)20. HPRD tissues were mapped to                   sets, the assumption of conditional independence required for a naive Bayes
                                                   the BRENDA Tissue Ontology21 (BTO) using direct matching where possible                  classifier is often not met, so we calculated and corrected for non-biological
                                                   and manual curation where direct matches were unavailable (Supplementary                 conditional dependency13.
                                                   Table 7). Tissues with fewer than ten directly annotated genes were pruned                   Each tissue model trained on the hierarchy-aware tissue-specific knowledge
                                                   as non-informative from 
```
### download
```
P values and tissue-                             can infer these features from large data compendia. We recently found
                                                   specific networks to identify disease-gene associations more                          that even samples measuring mixed cell lineages contain extractable
                                                   accurately than GWAS alone. Our webserver, GIANT, provides                            information related to lineage-specific expression 9. In addition to
                                                   an interface to human tissue networks through multi-gene                              tissue specificity, we10–13 and others14–17 have shown that hetero-
                                                   queries, network visualization, analysis tools including                              geneous genomic data contain functional information, for example,
                                                   NetWAS and downloadable networks. GIANT enables                                       of gene expression regulation by protein-DNA, protein-RNA,
                                                   systematic exploration of the landscape of interacting                                protein-protein and metabolite-protein interactions. Here we develop
npg




                                                   genes that shape specialized cellular functions across more                           and evaluate methods that simultaneously extract functional and
                                                   than a hundred human tissues and cell types.                                          tissue or cell-type signals to construct accurate maps of both where
                                                                                                                                         and how proteins act.
                                                   The precise actions of genes are frequently dependent on their tissue                    We build genome-scale functional maps of human tissues
                                                   context, and human diseases result from the disordered interplay of                   by integrating a collection of data sets covering thousands of
                                                                                                                                         experiments contained in more than 14,000 distinct publications.
                                                   1Department of Genetics, Geisel School of Medicine at Dartmouth, Hanover,
                                                                                                                                         To integrate these data, we automatically assess each data set for its
                                                   New Hampshire, USA. 2Dartmouth-Hitchcock Norris Cotton Cancer Center,
                                                   Lebanon, New Hampshire, USA. 3Institute for Quantitative Biomedical Sciences,
                                                                                                                                         relevance to each of 144 tissue- and cell lineage–specific functional
                                                   Dartmouth College, Hanover, New Hampshire, USA. 4Lewis-Sigler Institute for           contexts. The resulting functional maps provide a detailed portrait of
```
```
           research into disease mechanisms and therapy.                                                   ciation P values to receive NetWAS association scores. Visualizations
                                                                                                                                                   in the user-friendly, dynamic web interface are implemented using the
                                                   A dynamic, interactive interface for biomedical researchers                                     D3 library60, which enables use on any modern web browser without
                                                   To facilitate broad use of these networks by biomedical researchers,                            plugin installation. In addition to the interface, all of the underlying
                                                   we have developed GIANT—a dynamic, interactive web interface.                                   networks are provided for download, and the full list of input data sets
                                                   Researchers can query by individual genes or by gene sets of interest to                        and their sources is available through the webserver.
                                                   analyze tissue-specific gene functions and interactions. For example,
                                                   GIANT can provide tissue-specific functional maps and predictions                               DISCUSSION
                                                                                                                                                   Genes with tissue-specific expression and function have key roles in
                                                     a 0.8 Phenotypic: hypertension                                                                the physiological processes of complex organisms, and such genes
npg




                                                                                    disease gene prediction
                                                                            0.8                 0.8                     0.8                        are expected to underlie many human diseases2,3. Recent advances
                                                                          0.7          0.7          0.7                 0.7
                                                    OMIM (AUC)




                                                                          0.6          0.6          0.6                 0.6
                                                                                                                                                   Figure 5 Network reprioritization of hypertension GWAS identifies
                                                                                                                                                   hypertension-associated genes. Genes ranked using GWAS (gray) and
                                                                          0.5          0.5          0.5                 0.5                        genes reprioritized using NetWAS (brown) were assessed for
                                                                          0.4          0.4          0.4                 0.4
                                                                                                                                                   corresp
```
```
s. Am. J. Hum. Genet. 82, 949–958 (2008).
                                                       lymphoid enhancer binding factor 1 (LEF1). Development 126, 4441–4453                    64. Denny, J.C. et al. PheWAS: demonstrating the feasibility of a phenome-wide scan
                                                       (1999).                                                                                      to discover gene-disease associations. Bioinformatics 26, 1205–1210 (2010).
npg




                                                   576                                                                                                         VOLUME 47 | NUMBER 6 | JUNE 2015                 Nature Genetics
                                                   ONLINE METHODS                                                                           We defined ‘tissue categories’ from generic BRENDA terms, for example, nerv-
                                                   Data download and processing. We collected and integrated 987 genome-                    ous system, to categorize tissues into organ systems for evaluation and analysis.
                                                   scale data sets encompassing approximately 38,000 conditions from an                     For each tissue, we termed the set T as those genes directly annotated to that
                                                   estimated 14,000 publications including both expression and interaction                  tissue or any of its descendants in the ontology. We used tissue categories to
                                                   measurements. We downloaded interaction data from BioGRID65, IntAct66,                   define unrelated tissues (those not associated with the same category as the
                                                   MINT67 and MIPS68. BioGRID edges were discretized into five bins, labeled                tissue of interest). We defined T′ for each tissue as genes specifically annotated
                                                   0 to 4, where the bin number reflected the number of experiments supporting              to unrelated tissues.
                                                   the interaction. For the remaining databases, edges were discretized into the               Annotation of ubiquitously expressed genes. Genes ubiquitously expressed
                                                   presence or absence of an interaction.                                                   across tissues frequently carry out core biological processes and interact with
                                                      Predicting transcriptional regulation on the basis of DNA sequence is a               tissue-specific genes to perform specialized functions 78. We identified ubiq-
                                                   major challenge in understanding transcription at a systems level. To estimate           uitous genes from a multi-tissue RNA sequencing experiment79 and added
                                                   shared transcription factor regulation, binding motifs were downloaded from              ‘widely expressed’ genes from a multi–cell line mass spectroscopy experi-
                                                   JASPAR69. Genes were scored for the presence of transcription factor bind-               ment80, genes for proteins expressed in >75% of the tissues assayed in the
       
```
```
 NUMBER 6 | JUNE 2015                 Nature Genetics
                                                   ONLINE METHODS                                                                           We defined ‘tissue categories’ from generic BRENDA terms, for example, nerv-
                                                   Data download and processing. We collected and integrated 987 genome-                    ous system, to categorize tissues into organ systems for evaluation and analysis.
                                                   scale data sets encompassing approximately 38,000 conditions from an                     For each tissue, we termed the set T as those genes directly annotated to that
                                                   estimated 14,000 publications including both expression and interaction                  tissue or any of its descendants in the ontology. We used tissue categories to
                                                   measurements. We downloaded interaction data from BioGRID65, IntAct66,                   define unrelated tissues (those not associated with the same category as the
                                                   MINT67 and MIPS68. BioGRID edges were discretized into five bins, labeled                tissue of interest). We defined T′ for each tissue as genes specifically annotated
                                                   0 to 4, where the bin number reflected the number of experiments supporting              to unrelated tissues.
                                                   the interaction. For the remaining databases, edges were discretized into the               Annotation of ubiquitously expressed genes. Genes ubiquitously expressed
                                                   presence or absence of an interaction.                                                   across tissues frequently carry out core biological processes and interact with
                                                      Predicting transcriptional regulation on the basis of DNA sequence is a               tissue-specific genes to perform specialized functions 78. We identified ubiq-
                                                   major challenge in understanding transcription at a systems level. To estimate           uitous genes from a multi-tissue RNA sequencing experiment79 and added
                                                   shared transcription factor regulation, binding motifs were downloaded from              ‘widely expressed’ genes from a multi–cell line mass spectroscopy experi-
                                                   JASPAR69. Genes were scored for the presence of transcription factor bind-               ment80, genes for proteins expressed in >75% of the tissues assayed in the
                                                   ing sites using the MEME software suite70. FIMO71 was used to scan for each              human protein atlas81 and curated ‘ubiquitous genes’ from HPRD20. These
                                                   transcription factor profile within the 1-kb sequence upstream of each gene72.           8,475 ubiquitous genes (U) were considered expressed in all tissues and cell
                                                   Motif matches were treated as binary scores (present if P < 0.001). The final            types, in addition to the curated tissue-specific genes (T). Sets T and U were
                               
```
```
 to unrelated tissues.
                                                   the interaction. For the remaining databases, edges were discretized into the               Annotation of ubiquitously expressed genes. Genes ubiquitously expressed
                                                   presence or absence of an interaction.                                                   across tissues frequently carry out core biological processes and interact with
                                                      Predicting transcriptional regulation on the basis of DNA sequence is a               tissue-specific genes to perform specialized functions 78. We identified ubiq-
                                                   major challenge in understanding transcription at a systems level. To estimate           uitous genes from a multi-tissue RNA sequencing experiment79 and added
                                                   shared transcription factor regulation, binding motifs were downloaded from              ‘widely expressed’ genes from a multi–cell line mass spectroscopy experi-
                                                   JASPAR69. Genes were scored for the presence of transcription factor bind-               ment80, genes for proteins expressed in >75% of the tissues assayed in the
                                                   ing sites using the MEME software suite70. FIMO71 was used to scan for each              human protein atlas81 and curated ‘ubiquitous genes’ from HPRD20. These
                                                   transcription factor profile within the 1-kb sequence upstream of each gene72.           8,475 ubiquitous genes (U) were considered expressed in all tissues and cell
                                                   Motif matches were treated as binary scores (present if P < 0.001). The final            types, in addition to the curated tissue-specific genes (T). Sets T and U were
                                                   score for each gene pair was obtained by calculating the Pearson correlation             made disjoint by retaining only genes in T genes that were not in U.
                                                   between the motif association vectors for the genes.                                        Integration of tissue-specific and functional knowledge. We combined the
                                                      Chemical and genetic perturbation (c2:CGP) and microRNA target                        curated gene-to-tissue annotations with the tissue-naive functional gold stand-
                                                   (c3:MIR) profiles were downloaded from the Molecular Signatures Database                 ard to construct a hierarchical tissue-specific knowledgebase. We labeled each
                                                   (MSigDB73). Each gene pair’s score was the sum of shared profiles weighted               gene pair (positive or negative) in the functional relationship standard as spe-
                                                   by the specificity of each profile (1/len(genes)). The resulting scores were             cifically coexpressed in a tissue if both genes were tissue specific (T, T) or one
                                                   converted to z scores and discretized into bins ((−infinity, −1.5), [−1.5, −0.5),        was tissue specific and the other was ubiquitous (T, U). Interactions between
© 2015 Nature Am
```
```
s genes (U) were considered expressed in all tissues and cell
                                                   Motif matches were treated as binary scores (present if P < 0.001). The final            types, in addition to the curated tissue-specific genes (T). Sets T and U were
                                                   score for each gene pair was obtained by calculating the Pearson correlation             made disjoint by retaining only genes in T genes that were not in U.
                                                   between the motif association vectors for the genes.                                        Integration of tissue-specific and functional knowledge. We combined the
                                                      Chemical and genetic perturbation (c2:CGP) and microRNA target                        curated gene-to-tissue annotations with the tissue-naive functional gold stand-
                                                   (c3:MIR) profiles were downloaded from the Molecular Signatures Database                 ard to construct a hierarchical tissue-specific knowledgebase. We labeled each
                                                   (MSigDB73). Each gene pair’s score was the sum of shared profiles weighted               gene pair (positive or negative) in the functional relationship standard as spe-
                                                   by the specificity of each profile (1/len(genes)). The resulting scores were             cifically coexpressed in a tissue if both genes were tissue specific (T, T) or one
                                                   converted to z scores and discretized into bins ((−infinity, −1.5), [−1.5, −0.5),        was tissue specific and the other was ubiquitous (T, U). Interactions between
© 2015 Nature America, Inc. All rights reserved.




                                                   [–0.5, 0.5), [0.5, 1.5), [1.5, 2.5), [2.5, 3.5), [3.5, 4.5), [4.5, infinity)).           ubiquitous gene pairs were deemed not tissue specific and were ignored.
                                                      We downloaded all gene expression data sets from NCBI’s Gene Expression               After labeling specifically coexpressed gene pairs (edges) across all tissues,
                                                   Omnibus74 (GEO) and collapsed duplicate samples. GEO contains 980 human                  we considered four classes of edges—C1, C2, C3 and C4—to constitute each
                                                   data sets representing 20,868 conditions. Genes with more than 30% of values             tissue standard.
                                                   missing were removed, and remaining missing values were imputed using
                                                   ten neighbors75. Non-log-transformed data sets were log transformed.                       C1: positive functional edges between genes specifically coexpressed in
                                                   Expression measurements were summarized to Entrez76 identifiers, and                           the tissue [T–T and T–U].
                                                   duplicate identifiers were merged. The Pearson correlation was calculated for              C2: positive functional edges between a gene expressed in the tissue
                                                   each gene pair, normalized with Fisher’s z transform, mean subtracted and   
```
```
We labeled each
                                                   (MSigDB73). Each gene pair’s score was the sum of shared profiles weighted               gene pair (positive or negative) in the functional relationship standard as spe-
                                                   by the specificity of each profile (1/len(genes)). The resulting scores were             cifically coexpressed in a tissue if both genes were tissue specific (T, T) or one
                                                   converted to z scores and discretized into bins ((−infinity, −1.5), [−1.5, −0.5),        was tissue specific and the other was ubiquitous (T, U). Interactions between
© 2015 Nature America, Inc. All rights reserved.




                                                   [–0.5, 0.5), [0.5, 1.5), [1.5, 2.5), [2.5, 3.5), [3.5, 4.5), [4.5, infinity)).           ubiquitous gene pairs were deemed not tissue specific and were ignored.
                                                      We downloaded all gene expression data sets from NCBI’s Gene Expression               After labeling specifically coexpressed gene pairs (edges) across all tissues,
                                                   Omnibus74 (GEO) and collapsed duplicate samples. GEO contains 980 human                  we considered four classes of edges—C1, C2, C3 and C4—to constitute each
                                                   data sets representing 20,868 conditions. Genes with more than 30% of values             tissue standard.
                                                   missing were removed, and remaining missing values were imputed using
                                                   ten neighbors75. Non-log-transformed data sets were log transformed.                       C1: positive functional edges between genes specifically coexpressed in
                                                   Expression measurements were summarized to Entrez76 identifiers, and                           the tissue [T–T and T–U].
                                                   duplicate identifiers were merged. The Pearson correlation was calculated for              C2: positive functional edges between a gene expressed in the tissue
                                                   each gene pair, normalized with Fisher’s z transform, mean subtracted and                      and another specifically expressed in an unrelated tissue [T–T′ and
                                                   divided by the standard deviation. The resulting z scores were discretized                      U–T′].
                                                   into bins ((−infinity, −1.5), [–1.5, −0.5), [–0.5, 0.5), [0.5, 1.5), [1.5, 2.5), [2.5,     C3: negative functional edges between genes specifically coexpressed
                                                   3.5), [3.5, infinity)).                                                                        in the tissue [T–T and T–U].
                                                                                                                                              C4: negative functional edges between one gene expressed in the
                                                   Hierarchically aware knowledgebase construction via ontological pruning                        tissue and another specifically expressed in an unrelated tissue
                                                   with functional knowledge transf
```
### website
```
hropometric
                                                   additional GWAS: C-reactive protein levels (lnCRP)51, type 2 diabetes (T2D)88,              traits. PLoS Genet. 9, e1003500 (2013).
                                                                                                                                           90. Fritsche, L.G. et al. Seven new loci associated with age-related macular degeneration.
                                                   body mass index (BMI)89 and advanced age-related macular degeneration
                                                                                                                                               Nat. Genet. 45, 433–439 (2013).
                                                   (advanced AMD)90. Publicly available studies were obtained from their                   91. Mailman, M.D. et al. The NCBI dbGaP database of genotypes and phenotypes.
                                                   respective websites (BMI, advanced AMD) or the database of Genotypes and                    Nat. Genet. 39, 1181–1186 (2007).




                                                   doi:10.1038/ng.3259                                                                                                                                         Nature Genetics

```
### Supplementary Table 6
```
       with functional knowledge transfer.                                                            [T–T′ and U–T′].
                                                   Functional knowledge extraction. We constructed a tissue-naive functional
                                                   relationship gold standard from a set of 564 expert-selected GO biological               Among the four tissue classes, C1 represented tissue-specific functional
                                                   process terms and experimentally derived gene annotations (GO evidence                   relationships. To identify tissue-specific relationships, we constructed a
                                                   codes: EXP, IDA, IPI, IMP, IGI and IEP). Curators identified processes                   specific gold standard for each tissue by labeling edges in C1 as positives and
npg




                                                   testable through specific molecular experiments (Supplementary Table 6).                 edges in the other classes as negatives. Because C3 is defined on the basis of
                                                   Pairs of genes that were co-annotated to expert-selected terms after propaga-            tissue-expressed genes and C2 and C4 are defined on the basis of non-expressed
                                                   tion were treated as positive (i.e., functionally related) examples. Gene pairs          genes, the number of edges in these classes varied across tissues according
                                                   not co-annotated to any of these terms were considered as negative examples,             to how specific (cell type, tissue, organ or system), well studied (or easily
                                                   except in the following cases: (i) if two genes were annotated to two differ-            studied) and well curated (literature bias) they are. To construct comparable
                                                   ent GO terms with a significant number of shared genes (hypergeometric                   networks across tissues, we used a negative set composed of equal propor-
                                                   P value < 0.05) and (ii) if two genes were co-annotated to a set of ‘negative’           tions of edges from C2, C3 and C4. We limited all integrations to the set
                                                   GO terms that defined minimal relatedness77. Gene pairs that met either                  of 144 tissues (Supplementary Table 8) that contained at least ten C1 edges
                                                   condition were excluded from the set of negative examples and treated as                 between tissue-specific genes (T–T). This method incorporates the hierarchi-
                                                   neither related nor unrelated.                                                           cal relationships of tissues, allowing supervised methods to leverage these
                                                      Functional knowledge transfer. To increase the coverage of functional                 relationships.
                                                   interactions, we transferred experimentally confirmed mouse GO annota-                       Data integration. We constructed functional networks from genome-scale
                                                   tions to human functional analogs
```
### Supplementary Table 9
```
                                                                                                                        disease association. Mapping GO biological processes to tissues. To evaluate
                                                   Evaluation of tissue-specific functional relationships. We evaluated                tissue-specific functional rewiring in our networks, we needed associations
                                                   tissue-naive and tissue-specific functional networks using fivefold cross-          between tissues and tissue-specific processes. We used text matching followed
                                                   validation. The 6,062 genes represented in the tissue-specific knowledge-           by manual curation to map biological process (BP) terms in GO to tissue terms
                                                   base were randomly partitioned into 5 sets. For each cross-validation run,          in the BRENDA Tissue ontology (Supplementary Table 9).
                                                   gene pairs where neither gene was present in the holdout interval were used             Network connectivity of tissue-specific processes. For each tissue, we
                                                   for training. Any gene pair where both genes were present in the holdout            constructed a tissue-minus-naive network by subtracting edge probabilities
                                                   was used for evaluation of the AUC. The estimated performance of each of            of the naive network from those of the tissue network. Negative weights were
                                                   the 144 functional networks was summarized as the median AUC of the five            set to zero. In this subtracted network, positive scores corresponded to edges
                                                   cross-validation runs (Supplementary Table 8).                                      with a tissue network interaction probability greater than the naive network
                                                      Mapping data sets to tissues. We mapped data sets to tissues to compare with     probability. We expected relevant tissue-specific processes to be more con-
                                                   an integration of only tissue-specific data. On the basis of previous work83 that   nected in the tissue network than the naive network and over processes that
                                                   annotated samples from biological text, we extracted the title and description      are not. For instance, for T lymphocytes, ‘T cell receptor signaling pathway’
                                                   for each GDS data set and annotated each using MetaMap84. This resulted in a        is a relevant process, whereas ‘neuron projection development’ is not. Within
                                                   mapping of GDS data sets to Unified Medical Language System (UMLS) terms.           each subtracted tissue network, we ranked all tissue-specific processes by
                                                   We applied the same process for the title and description of each BRENDA            their edge density in the network and evaluated the extent to which relevant
                                                   tissue and merged the two mappings by shared UMLS terms.                            processes (positives) were ranked 
```
## GraphSAGE.txt

### URLs
- `http://snap.stanford.edu/graphsage/`
- `http://snap.stanford.edu/graphsage/.`
- `https://archive.`
- `https://archive.org/details/`
- `https://github.com/tensorflow/models/blob/master/tutorials/embedding/`
### Gene Ontology
```
ains compared


                                                        7
to the GCN approach. For example, the unsupervised variant GraphSAGE-pool outperforms the
concatenation of the DeepWalk embeddings and the raw features by 13.8% on the citation data
and 29.1% on the Reddit data, while the supervised version provides a gain of 19.7% and 37.2%,
respectively. Interestingly, the LSTM based aggregator shows strong performance, despite the fact
that it is designed for sequential data and not unordered sets. Lastly, we see that the performance of
unsupervised GraphSAGE is reasonably competitive with the fully supervised version, indicating
that our framework can achieve strong performance without task-specific fine-tuning.

4.2       Generalizing across graphs: Protein-protein interactions

We now consider the task of generalizing across graphs, which requires learning about node roles
rather than community structure. We classify protein roles—in terms of their cellular functions from
gene ontology—in various protein-protein interaction (PPI) graphs, with each graph corresponding
to a different human tissue [41]. We use positional gene sets, motif gene sets and immunological
signatures as features and gene ontology sets as labels (121 in total), collected from the Molecular
Signatures Database [34]. The average graph contains 2373 nodes, with an average degree of 28.8.
We train all algorithms on 20 graphs and then average prediction F1 scores on two test graphs (with
two other graphs used for validation).
The final two columns of Table 1 summarize the accuracies of the various approaches on this
data. Again we see that GraphSAGE significantly outperforms the baseline approaches, with the
LSTM- and pooling-based aggregators providing substantial gains over the mean- and GCN-based
aggregators.6

4.3       Runtime and parameter sensitivity

Figure 2.A summarizes the training and test runtimes for the different approaches. The training time
for the methods are comparable (with GraphSAGE-LSTM being the slowest). However, the need to
sample new random walks and run new rounds of SGD to embed unseen nodes makes DeepWalk
100-500× slower at test time.
For the GraphSAGE variants, we found that setting K = 2 provided a consistent boost in accuracy of
around 10-15%, on average, compared to K = 1; however, increasing K beyond 2 gave marginal
returns in performance (0-5%) while increasing the runtime by a prohibitively large factor of 10-100×,
depending on the neighborhood sample size. We also found diminishing returns for sampling
large neighborhoods (Figure 2.B). Thus, despite the higher variance induced by sub-sampling
neighborhoods, GraphSAGE is still able to maintain strong predictive accuracy, while significantly
improving the runtime.

4.4       Summary comparison between the different aggregator architectures

Overall, we found that the LSTM- and pool-based aggregators performed the best, in terms of both
average performance and number of experimental settings where they were the top-performing
method (Table 1). To give more quantitative insight into these trends, we consider each of the
six different experimental settings (i.e., (3 datasets) × (unsupervised vs. supervised)) as trials and
consider what performance trends are likely to generalize. In particular, we use the non-parametric
Wilcoxon Signed-Rank Test [33] to quantify the differences between the different aggregators across
trials, reporting the T -statistic and p-value where applicable. Note t
```
```
tures by 13.8% on the citation data
and 29.1% on the Reddit data, while the supervised version provides a gain of 19.7% and 37.2%,
respectively. Interestingly, the LSTM based aggregator shows strong performance, despite the fact
that it is designed for sequential data and not unordered sets. Lastly, we see that the performance of
unsupervised GraphSAGE is reasonably competitive with the fully supervised version, indicating
that our framework can achieve strong performance without task-specific fine-tuning.

4.2       Generalizing across graphs: Protein-protein interactions

We now consider the task of generalizing across graphs, which requires learning about node roles
rather than community structure. We classify protein roles—in terms of their cellular functions from
gene ontology—in various protein-protein interaction (PPI) graphs, with each graph corresponding
to a different human tissue [41]. We use positional gene sets, motif gene sets and immunological
signatures as features and gene ontology sets as labels (121 in total), collected from the Molecular
Signatures Database [34]. The average graph contains 2373 nodes, with an average degree of 28.8.
We train all algorithms on 20 graphs and then average prediction F1 scores on two test graphs (with
two other graphs used for validation).
The final two columns of Table 1 summarize the accuracies of the various approaches on this
data. Again we see that GraphSAGE significantly outperforms the baseline approaches, with the
LSTM- and pooling-based aggregators providing substantial gains over the mean- and GCN-based
aggregators.6

4.3       Runtime and parameter sensitivity

Figure 2.A summarizes the training and test runtimes for the different approaches. The training time
for the methods are comparable (with GraphSAGE-LSTM being the slowest). However, the need to
sample new random walks and run new rounds of SGD to embed unseen nodes makes DeepWalk
100-500× slower at test time.
For the GraphSAGE variants, we found that setting K = 2 provided a consistent boost in accuracy of
around 10-15%, on average, compared to K = 1; however, increasing K beyond 2 gave marginal
returns in performance (0-5%) while increasing the runtime by a prohibitively large factor of 10-100×,
depending on the neighborhood sample size. We also found diminishing returns for sampling
large neighborhoods (Figure 2.B). Thus, despite the higher variance induced by sub-sampling
neighborhoods, GraphSAGE is still able to maintain strong predictive accuracy, while significantly
improving the runtime.

4.4       Summary comparison between the different aggregator architectures

Overall, we found that the LSTM- and pool-based aggregators performed the best, in terms of both
average performance and number of experimental settings where they were the top-performing
method (Table 1). To give more quantitative insight into these trends, we consider each of the
six different experimental settings (i.e., (3 datasets) × (unsupervised vs. supervised)) as trials and
consider what performance trends are likely to generalize. In particular, we use the non-parametric
Wilcoxon Signed-Rank Test [33] to quantify the differences between the different aggregators across
trials, reporting the T -statistic and p-value where applicable. Note that this method is rank-based and
essentially tests whether we would expect one particular approach to outperform another in a new
experimental setting. Given our small sample size of only 6 different settings, this signi
```
### website
```
e edges in all graphs before feeding them into the GraphSAGE algorithm. In particular,
we subsample edges so that no node has degree larger than 128. Since we only sample at most 25
neighbors per node, this is a reasonable tradeoff. This downsampling allows us to store neighborhood
information as dense adjacency lists, which drastically improves computational efficiency. For the
Reddit data we also downsampled the edges of the original graph as a pre-processing step, since the
    7
      Note that these values differ from our previous reported pre-print values because they are corrected to account
for an extraneous normalization by the batch size. We thank Ben Johnson for pointing out this discrepancy.
    8
      https://github.com/tensorflow/models/blob/master/tutorials/embedding/
word2vec.py


                                                         14
original graph is extremely dense. All experiments are on the downsampled version, but we release
the full version on the project website for reference.


D    Alignment Issues and Orthogonal Invariance for DeepWalk and Related
     Approaches

DeepWalk [28], node2vec [11], and other recent successful node embedding approaches employ
objective functions of the form:
                                   X                X
                                 α   f (z>
                                         i zj ) + β    g(z>
                                                          i zj )                       (4)
                                   i,j∈A               i,j∈B

where f , g are smooth, continuous functions, zi are the node representations that are being directly
optimized (i.e., via embedding look-ups), and A, B are sets of pairs of nodes. Note that in many cases,
in the actual code implementations used by the authors of these approaches, nodes are associated
with two unique embedding vectors and the arguments to the dot products in f and g are drawn for
distinct embedding look-ups (e.g., [11, 28]); however, this does not fundamentally alter the learning
algorithm. The majority of approaches also normalize the learned embeddings to unit length, so we
assume this post-processing as well.
By connection to word embedding approaches and the arguments of [20], these approaches can
also be viewed as stochastic, implicit matrix factorizations where we are trying to learn a matrix
Z ∈ R|V|×d such that
                                           ZZ> ≈ M,                                            (5)
where M is some matrix containing random walk statistics.
An important consequence of this structure is that the embeddings can be rotated by an arbitrary
orthogonal matrix, without impacting the objective:

                                           ZQ> QZ> = ZZ> ,                                         (6)
where Q ∈ Rd×d is any orthogonal matrix. Since the embeddings are otherwise unconstrained and
the only error signal comes from the orthogonally-invariant objective (4), the entire embedding space
is free to arbitrarily rotate during training.
Two clear consequences of this are:

      1. Suppose we run an embedding approach based on (4) on two separate graphs A and B
         using the same output dimension. Without some explicit penalty enforcing alignment, the
         learned embeddings spaces for the two graphs will be arbitrarily rotated with respect to each
         other after training. Thus, for any node classification method that is trained on individual
         embeddings from gra
```